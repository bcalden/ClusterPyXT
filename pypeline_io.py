import os
import errno
import shutil
import gzip
import glob
import re
import pycrates as pc
import sys
import csv

class Colors:
    BLACK = "\u001b[30m"
    RED = "\u001b[31m"
    GREEN = "\u001b[32m"
    YELLOW = "\u001b[33m"
    BLUE = "\u001b[34m"
    MAGENTA = "\u001b[35m"
    CYAN = "\u001b[36m"
    WHITE = "\u001b[37m"
    RESET = "\033[0;0m"
    BOLD = "\033[;1m"

class Ansi:

    # clearing screen codes
    CLEAR_UNTIL_END_OF_SCREEN = "\u001b[0J"
    CLEAR_UNTIL_START_OF_SCREEN = "\u001b[1J"
    CLEAR_SCREEN = "\u001b[2J"

    # clearing line codes
    CLEAR_UNTIL_END_OF_LINE = "\u001b[0K"
    CLEAR_UNTIL_START_OF_LINE = "\u001b[1K"
    CLEAR_LINE = "\u001b[2K"

def flush():
    sys.stdout.flush()

def write(string):
    sys.stdout.write(string)

def clear_line():
    sys.stdout.write(Ansi.CLEAR_LINE)
    move_cursor_left(1000)
    return

def reset_cursor():
    move_cursor_left(1000)
    move_cursor_up(1000)

def move_cursor_up(number_spaces):
    sys.stdout.write(u"\u001b[{n:d}A".format(n=number_spaces))
    return


def move_cursor_down(number_spaces):
    sys.stdout.write(u"\u001b[{n:d}B".format(n=number_spaces))
    return


def move_cursor_right(number_spaces):
    sys.stdout.write(u"\u001b[{n:d}C".format(n=number_spaces))
    return


def move_cursor_left(number_spaces):
    sys.stdout.write(u"\u001b[{n:d}D".format(n=number_spaces))
    return


def clear_screen():
    sys.stdout.write(Ansi.CLEAR_SCREEN)
    reset_cursor()
    sys.stdout.flush()


def get_pixel_values(filename):
    image = pc.read_file(filename)
    return pc.copy_piximgvals(image) # returns a numpy array


def make_directory(directory):
    try:
        os.makedirs(get_path(directory))
    except OSError as err:
        if err.errno == errno.EEXIST:
            print("Directory {} already exists, skipping.".format(directory))
            pass
        elif err.errno in [errno.EACCES, errno.EROFS]:
            print("Unable to create {}, check user permissions allow for write access to {} and its"
                  "parent directories.".format(directory, directory))
            raise
        elif err.errno == errno.ENOSPC:
            print("No space left to create {}.".format(directory))
            raise
        else:
            raise


def make_initial_directories(cluster_obj):
    # Questioning if this function needs to exist...
    print("Making directories")
    directories_to_make = ['combined', 'wvt', 'acb', 'main_output']
    cluster_dir = cluster_obj.directory
    full_paths = [get_path("{}/{}".format(cluster_dir, directory)) for directory in directories_to_make]
    full_paths += [get_path("{}/{}/analysis/".format(cluster_dir, obsid)) for obsid in cluster_obj.observation_ids]
    for directory in full_paths:
        make_directory(directory)
    return


def set_working_directory(directory):
    try:
        os.chdir(directory)
    except NotADirectoryError as err:
        print(err)
        print("Error: Not a directory")
        raise


def get_user_input(prompt, data_description=""):
    user_input_good = False
    user_input = None
    while not user_input_good:
        user_input = input(prompt)
        if data_description != "":
            y_no_flag = input("Using {} for {}. Is this correct? [y/n]: ".format(user_input, data_description))
        else:
            y_no_flag = input("You entered {}. Is this correct? [y/n]:".format(user_input))

        if y_no_flag.lower() in ['y', 'yes']:
            user_input_good = True

    return user_input

# To add: A I/O function for getting and verifying an input directory/file

def check_yes_no(prompt):
    y_no_flag = False
    valid_input = False
    while not valid_input:
        # user_input = input("{}\n Is this correct? (y/n): ".format(prompt))
        user_input = input(prompt)
        if user_input.lower() in ['y', 'yes']:
            y_no_flag = True
            valid_input = True
        elif user_input.lower() in ['n', 'no']:
            y_no_flag = False
            valid_input = True

    return y_no_flag


def get_path(path):
    return os.path.normpath(path)


def write_contents_to_file(contents, filename, binary=True):
    try:
        file_attributes = 'wb' if binary else 'w'
        with open(filename, file_attributes) as f:
            f.write(contents)
    except FileExistsError:
        print("File already exists and I can't over write!")
        raise
    except FileNotFoundError:
        print("File not found!")
        raise
    return


def read_contents_of_file(filename):
    try:
        with open(filename) as f:
            data = f.read()
    except FileNotFoundError:
        print("File %s Not Found!" % filename)
        raise

    return data


def copy(src, dst):
    shutil.copy2(src, dst)


def copytree(src, dst):  # , symlinks=False, ignore=None):
    # for item in os.listdir(src):
    #     s = os.path.join(src, item)
    #     d = os.path.join(dst, item)
    #     if os.path.isdir(s):
    #         shutil.copytree(s, d, symlinks, ignore)
    #     else:
    #         shutil.copy2(s, d)

    # http://stackoverflow.com/questions/7419665/python-move-and-overwrite-files-and-folders

    for src_dir, dirs, files in os.walk(src):
        dst_dir = src_dir.replace(src, dst, 1)
        if not os.path.exists(dst_dir):
            os.makedirs(dst_dir)
        for file_ in files:
            src_file = os.path.join(src_dir, file_)
            dst_file = os.path.join(dst_dir, file_)
            if os.path.exists(dst_file):
                os.remove(dst_file)
            shutil.copy(src_file, dst_dir)

    return


def gz_unzip(source): #, destination=""):
    input_file = gzip.open(source, 'rb')
    output_file_name = source[:-3]
    output_file = open(output_file_name, 'wb')
    output_file.write(input_file.read())
    input_file.close()
    output_file.close()

    return


def get_filename_matching(pattern_str):
    return glob.glob(pattern_str)


def get_date_from_filename(filename):
    # acisD2000-01-29gain_ctiN0006.fits

    regex = r"(19|20)\d\d[- /.](0[1-9]|1[012])[- /.](0[1-9]|[12][0-9]|3[01])"

    return get_regex_result_from_string(regex, filename)


def get_regex_result_from_string(regex, string):
    match = re.search(regex, string)

    return match.group()

def get_version_from_filename(filename):
    # acisD2000-01-29gain_ctiN0006.fits

    regex = r"N\d{0,}."

    return get_regex_result_from_string(regex, filename)


def grep(filename, grep):
    with open(filename) as f:
        data = f.readlines()
    data = [x.strip() for x in data]

    found = []

    for line in data:
        if -1 != line.find("_pnt"):
            found.append(line)

    return found


def change_extension(filename, new_extension):
    split_filename = filename.split('.')
    split_filename[-1] = new_extension
    return '.'.join(split_filename)


def get_sources_filename():
    prompt = "Please enter the full path for the sources.reg file: "
    filename = get_user_input(prompt)

    return filename


def file_exists(filename):
    return os.path.isfile(filename)


def file_size(filename):
    return os.path.getsize(filename)


def delete(filename):
    os.remove(filename)


def delete_if_exists(filename):
    if file_exists(filename):
        print("Deleting {filename}.".format(filename=filename))
        delete(filename)
    return


def read_line_number(filename, line_number):
    if line_number < 1:
        return None
    for current_line_number, line in enumerate(open(filename, 'rU')):
        if current_line_number == line_number - 1:
            return line[:-1] if line[-1:] == '\n' else line
    return None


def append_to_file(filename, data):
    with open(filename, 'a') as f:
        f.write(data)


def num_lines_in(filename):
    with open(filename) as f:
        for i, l in enumerate(f):
            pass
        return i+1


def get_cluster_info_from_csv(csv_file):
    clusters = []

    with open(csv_file, 'r') as f:
        reader = csv.DictReader(f)
        for row in reader:
            name = row['name']
            obsids = [str(x) for x in row['obsids'].split(' ')]
            abundance = float(row['abundance'])
            nH = float(row['nH'])
            z = float(row['z'])
            clusters.append({'name': name,
                             'obsids': obsids,
                             'abundance': abundance,
                             'hydrogen_column_density': nH,
                             'redshift': z})
    return clusters

def make_initial_data_dir(directory):
    if not os.path.isdir(directory):
        make_directory(directory)