import config
import glob
import pypeline_io as io
from enum import Enum
import curses
import sys
import os
import cluster
import ciao


class MenuItemType(Enum):
    TITLE = 0
    SUBTITLE = 1
    SUBMENU = 2
    COMMAND = 3
    FUNCTION = 4
    CLUSTER = 5

    SYSTEMEXIT = -1


class MenuOption:
    def __init__(self, title="", type="", value="", enabled=True):
        self.title = title
        self.type = type
        self.value = value
        self.enabled = enabled


class CursesMenu(object):

    INIT = {'type' : 'init'}

    def __init__(self, menu_options):
        self.screen = curses.initscr()
        self.menu_options = menu_options
        self.selected_option = 0
        self._previously_selected_option = None
        self.running = True

        #init curses and curses input
        curses.noecho()
        curses.cbreak()
        curses.start_color()
        curses.curs_set(0) #Hide cursor
        self.screen.keypad(1)

        #set up color pair for highlighted option
        curses.init_pair(1, curses.COLOR_BLACK, curses.COLOR_WHITE)
        self.hilite_color = curses.color_pair(1)
        self.normal_color = curses.A_NORMAL

        self.title_line_buffer = \
            len(self.menu_options['title']) + 3 if isinstance(self.menu_options['title'], list) else 1

    def prompt_selection(self, parent=None):
        if parent is None:
            lastoption = "Exit"
        else:
            lastoption = "Return to previous menu ({})".format(parent['title'])

        option_count = len(self.menu_options['options'])

        input_key = None

        ENTER_KEY = ord('\n')
        while input_key != ENTER_KEY:
            if self.selected_option != self._previously_selected_option:
                self._previously_selected_option = self.selected_option

            self.screen.border(0)
            self._draw_title()
            for option in range(option_count):
                if self.selected_option == option:
                    self._draw_option(option, self.hilite_color)
                else:
                    self._draw_option(option, self.normal_color)

            if self.selected_option == option_count:
                self.screen.addstr(self.title_line_buffer + option_count+2, 4, "{:2} - {}".format(option_count+1,
                    lastoption), self.hilite_color)
            else:
                self.screen.addstr(self.title_line_buffer + option_count+2, 4, "{:2} - {}".format(option_count+1,
                    lastoption), self.normal_color)

            max_y, max_x = self.screen.getmaxyx()
            if input_key is not None:
                self.screen.addstr(max_y-3, max_x - 5, "{:3}".format(self.selected_option))
            self.screen.refresh()


            input_key = self.screen.getch()
            down_keys = [curses.KEY_DOWN, ord('j')]
            up_keys = [curses.KEY_UP, ord('k')]
            exit_keys = [ord('q')]

            if input_key in down_keys:
                if self.selected_option < option_count:
                    self.selected_option += 1
                else:
                    self.selected_option = 0

            if input_key in up_keys:
                if self.selected_option > 0:
                    self.selected_option -= 1
                else:
                    self.selected_option = option_count

            if input_key in exit_keys:
                self.selected_option = option_count #auto select exit and return
                break

        return self.selected_option

    def _draw_option(self, option_number, style):
        self.screen.addstr(self.title_line_buffer+2 + option_number,
                           4,
                           "{:2} - {}".format(option_number+1, self.menu_options['options'][option_number]['title']),
                           style)

    def _draw_title(self):
        title_str = self.menu_options['title']

        for i, line in enumerate(title_str):
            self.screen.addstr(2+i, 2, line, curses.A_NORMAL)
        self.screen.addstr(self.title_line_buffer, 2, self.menu_options['subtitle'], curses.A_BOLD)

    def display(self):
        selected_option = self.prompt_selection()
        i, _ = self.screen.getmaxyx()
        #curses.endwin()
        os.system('clear')
        if selected_option < len(self.menu_options['options']):
            selected_opt = self.menu_options['options'][selected_option]
            return selected_opt
        else:
            self.running = False
            return {'title': 'Exit', 'type': 'exitmenu'}


def get_cluster_name_from_config_file(config_file):
    return os.path.basename(os.path.split(config_file)[0])


def get_cluster_configs(data_dir=config.data_directory()):
    dirs = os.listdir(data_dir)
    configuration_files = []
    cluster_names = []
    for directory in dirs:
        #print("{}/{}".format(data_dir,directory))
        config_file = glob.glob("{data_dir}/{directory}/*_pypeline_config.ini".format(
            data_dir=data_dir,
            directory=directory
        ))
        if config_file:
            configuration_files.append(config_file[0])
            cluster_names.append(get_cluster_name_from_config_file(config_file[0]))
    return list(zip(cluster_names, configuration_files))


# def last_cluster_worked_on():
#     current_cluster_file = config.current_cluster_file()
#     current_cluster_name = config.current_cluster_name()
#
#     return current_cluster_name, current_cluster_file

def continue_cluster():
    current_cluster = config.current_cluster()
    if current_cluster is not None:
        ciao.start_from_last(current_cluster)
    else:
        print("No current cluster / Unable to find current cluster."
              "Check the .ini files in the cluster and pypeline directories")


def change_current_cluster(cluster_config):
    new_cluster = cluster.read_cluster_data(cluster_config)
    config.set_current_cluster(new_cluster)


def initialize_cluster():
    new_cluster = cluster.ClusterObj()
    new_cluster.initialize_cluster()
    config.set_current_cluster(new_cluster)


def display_menu(menu):
    selected_action = menu.display()
    if selected_action['type'] != 'exitmenu':
        if selected_action['type'] == MenuItemType.COMMAND:
            os.system(selected_action['command'])
        elif selected_action['type'] == MenuItemType.FUNCTION:
            curses.nocbreak()
            menu.screen.keypad(True)
            curses.echo()
            curses.endwin()
            selected_action['function']()
        elif selected_action['type'] == MenuItemType.SUBMENU:
            display_menu(selected_action['menu'])
        elif selected_action['type'] == MenuItemType.CLUSTER:
            change_current_cluster(selected_action['configuration'])
    else:
        pass
#        curses.reset_shell_mode()


def main_menu_builder():
    current_cluster = config.current_cluster()
    options = []

    subtitle_string = "Current working cluster: {cluster}"
    if current_cluster is not None:
        subtitle = subtitle_string.format(cluster=current_cluster.name)

        continue_cluster_string = 'Continue {cluster}'.format(
            cluster=current_cluster.name
        )

        options.append({'title': continue_cluster_string,
                        'type': MenuItemType.FUNCTION,
                        'function': continue_cluster
                        })
    else:
        subtitle = subtitle_string.format(cluster="None.")

    menu = {'title': fancy_title_3(),
            'type': 'menu',
            'subtitle': subtitle}

    options.append({'title': 'Initialize a new cluster',
                    'type': MenuItemType.FUNCTION,
                    'function': initialize_cluster})
    options.append({'title': 'Change current cluster',
                    'type': MenuItemType.SUBMENU,
                    'menu': change_cluster_menu()})

    menu['options'] = options

    return CursesMenu(menu)


def change_cluster_menu():
    menu = {'title': fancy_title(),
            'type': MenuItemType.SUBMENU,
            'subtitle': 'Choose a cluster'}
    options = []
    current_clusters = get_cluster_configs()
    for cluster in current_clusters:
        options.append({
            'title': cluster[0],
            'type': MenuItemType.CLUSTER,
            'configuration': cluster[1]
        })
    menu['options'] = options

    return CursesMenu(menu)

def fancy_title():
    title = ["  __  __                       ____                   _ _",
             "  \ \/ /     _ __ __ _ _   _  |  _ \ _   _ _ __   ___| (_)_ __   ___ ",
             "   \  /_____| '__/ _` | | | | | |_) | | | | '_ \ / _ \ | | '_ \ / _ \\",
             "   /  \_____| | | (_| | |_| | |  __/| |_| | |_) |  __/ | | | | |  __/",
             "  /_/\_\    |_|  \__,_|\__, | |_|    \__, | .__/ \___|_|_|_| |_|\___|",
             "                       |___/         |___/|_|                        "]

    return title

def fancy_title_2():
    title = [" ██████╗  █████╗ ██╗      █████╗ ██╗  ██╗██╗   ██╗     ██████╗██╗     ██╗   ██╗███████╗████████╗███████╗██████╗ ",
             "██╔════╝ ██╔══██╗██║     ██╔══██╗╚██╗██╔╝╚██╗ ██╔╝    ██╔════╝██║     ██║   ██║██╔════╝╚══██╔══╝██╔════╝██╔══██╗",
             "██║  ███╗███████║██║     ███████║ ╚███╔╝  ╚████╔╝     ██║     ██║     ██║   ██║███████╗   ██║   █████╗  ██████╔╝",
             "██║   ██║██╔══██║██║     ██╔══██║ ██╔██╗   ╚██╔╝      ██║     ██║     ██║   ██║╚════██║   ██║   ██╔══╝  ██╔══██╗",
             "╚██████╔╝██║  ██║███████╗██║  ██║██╔╝ ██╗   ██║       ╚██████╗███████╗╚██████╔╝███████║   ██║   ███████╗██║  ██║",
             " ╚═════╝ ╚═╝  ╚═╝╚══════╝╚═╝  ╚═╝╚═╝  ╚═╝   ╚═╝        ╚═════╝╚══════╝ ╚═════╝ ╚══════╝   ╚═╝   ╚══════╝╚═╝  ╚═╝",
             "                                                                                                                ",
             "████████╗███╗   ███╗ █████╗ ██████╗         ██████╗ ██╗   ██╗██████╗ ███████╗██╗     ██╗███╗   ██╗███████╗      ",
             "╚══██╔══╝████╗ ████║██╔══██╗██╔══██╗        ██╔══██╗╚██╗ ██╔╝██╔══██╗██╔════╝██║     ██║████╗  ██║██╔════╝      ",
             "   ██║   ██╔████╔██║███████║██████╔╝        ██████╔╝ ╚████╔╝ ██████╔╝█████╗  ██║     ██║██╔██╗ ██║█████╗        ",
             "   ██║   ██║╚██╔╝██║██╔══██║██╔═══╝         ██╔═══╝   ╚██╔╝  ██╔═══╝ ██╔══╝  ██║     ██║██║╚██╗██║██╔══╝        ",
             "   ██║   ██║ ╚═╝ ██║██║  ██║██║             ██║        ██║   ██║     ███████╗███████╗██║██║ ╚████║███████╗      ",
             "   ╚═╝   ╚═╝     ╚═╝╚═╝  ╚═╝╚═╝             ╚═╝        ╚═╝   ╚═╝     ╚══════╝╚══════╝╚═╝╚═╝  ╚═══╝╚══════╝      "]
    return title

def fancy_title_3():
    title = [" _____ _           _           ______     __   _______ ",
             "/  __ \ |         | |          | ___ \    \ \ / /_   _|",
             "| /  \/ |_   _ ___| |_ ___ _ __| |_/ /   _ \ V /  | |  ",
             "| |   | | | | / __| __/ _ \ '__|  __/ | | |/   \  | |  ",
             "| \__/\ | |_| \__ \ ||  __/ |  | |  | |_| / /^\ \ | |  ",
             " \____/_|\__,_|___/\__\___|_|  \_|   \__, \/   \/ \_/  ",
             "                                      __/ |            ",
             "                                     |___/             ",
             "The Galaxy Cluster Pypeline for X-ray Temperature Maps "]

    return title


def make_menu():
    main_menu = main_menu_builder()
    display_menu(main_menu)
    io.move_cursor_left(1000)
    print()
    print("Be sure to note any instructions before ending.")
    io.move_cursor_left(1000)
    print("Press enter to end... (screen will be erased)")
    main_menu.screen.getch()
    curses.nocbreak()
    main_menu.screen.keypad(False)
    curses.echo()
    curses.endwin()


if __name__=="__main__":
    make_menu()