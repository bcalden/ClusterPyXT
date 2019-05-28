import os
import cluster
import pypeline_io as io
import configparser

CONFIG_FILENAME = 'pypeline_config.ini'

###
# rewrite to make use of ciao environment variables for CALDB and CIAO dir

def get_user_input(prompt):
    user_input_good = False
    user_input = None
    while not user_input_good:
        user_input = input(prompt)
        y_no_flag = input("You entered {}. Is this correct? (y/n):".format(user_input))

        if y_no_flag.lower() in ['y', 'yes']:
            user_input_good = True

    return user_input


def set_current_cluster(cluster):
    update_system_configuration(current_cluster_name=cluster.name,
                                cluster_config_file=cluster.configuration_filename)


def initialize_system_configuration(filename):
    config = configparser.ConfigParser()

    data_directory = get_user_input("What directory should be used for data storage? (Please enter the full path): ")
    ciao_directory = get_user_input("What directory is CIAO located in? (Please enter the full path,"
                                    " for exampe: /somedir/otherdir/ciao/ciao-4.9): ")

    io.make_initial_data_dir(data_directory)

    config_dict = {'data_directory': os.path.normpath(data_directory),
                   'ciao_directory': os.path.normpath(ciao_directory),
                   'current_cluster_name': '',
                   'cluster_config_file': ''}

    config['pypeline_config'] = config_dict

    with open(filename, 'w') as configfile:
        config.write(configfile)

    return


def update_system_configuration(**kwargs):
    current_config = read_system_configuration(CONFIG_FILENAME)

    for key in ('data_directory', 'ciao_directory', 'current_cluster_name', 'cluster_config_file'):
        if key in kwargs:
            current_config[key] = kwargs[key]

    new_config = configparser.ConfigParser()

    new_config['pypeline_config'] = current_config

    with open(CONFIG_FILENAME, 'w') as configfile:
        new_config.write(configfile)

    return


def read_system_configuration(filename):
    config = configparser.ConfigParser()

    config.read(filename)
    try:
        config_dict = dict(config['pypeline_config'])
    except KeyError as err:
        print(err)
        raise
    else:
        if ("" != config_dict['current_cluster_name']) or ("" == config_dict['cluster_config_file']):
            cluster_folder_exists = os.path.exists("{}/{}".format(config_dict['data_directory'],
                                                              config_dict['current_cluster_name']))
            cluster_config_exists = os.path.exists(config_dict['cluster_config_file'])

            if not (cluster_config_exists and cluster_folder_exists):
                config_dict['current_cluster_name'] = ""
                config_dict['cluster_config_file'] = ""
                print("No cluster currently being processed.")

    return config_dict


def system_configuration():
    return read_system_configuration(CONFIG_FILENAME)


def data_directory():
    return system_configuration()['data_directory']


def ciao_directory():
    return system_configuration()['ciao_directory']


def current_cluster_name():
    current_cluster = system_configuration()['current_cluster_name']

    if current_cluster in [None, '']:
        current_cluster = None

    return current_cluster


def current_cluster_file():
    current_file = system_configuration()['cluster_config_file']

    if current_file in [None, '']:
        current_file = None

    return current_file


def current_cluster():
    if current_cluster_file() is not None:
        return cluster.read_cluster_data(current_cluster_file())
    else:
        return None

def initialize_pypeline():
    system_config = None

    while system_config is None:
        try:
            system_config = read_system_configuration(CONFIG_FILENAME)
        except KeyError as err:
                initialize_system_configuration(CONFIG_FILENAME)
        except:
            print("Unknown error.")
            raise

    #add_ciao_to_path(system_config['ciao_directory'])
