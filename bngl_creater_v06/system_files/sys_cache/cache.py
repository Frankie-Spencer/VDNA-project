import os
# from sys_cache.search_chrome_browser import get_bowser_location


def read_s_loc():
    with open('system_files/sys_cache/_source_loc.cache', 'r') as cf:
        li = cf.readline()
        cf.close()
        return li


def write_s_loc(l):
    with open('system_files/sys_cache/_source_loc.cache', 'w') as cf:
        li = os.path.dirname(l)
        cf.writelines(li)
        cf.close()


def read_d_loc():
    with open('system_files/sys_cache/_des_loc.cache', 'r') as cf:
        li = cf.readline()
        cf.close()
        return li


def write_d_loc(l):
    with open('system_files/sys_cache/_des_loc.cache', 'w') as cf:
        cf.writelines(l)
        cf.close()


def write_browser_loc(l):
    with open('system_files/sys_cache/_browser_loc.cache', 'w') as cf:
        # print(l, 'write call')
        cf.writelines(l)
        cf.close()


def read_browser_loc():
    with open('system_files/sys_cache/_browser_loc.cache', 'r') as cf:
        li = cf.readline()
        cf.close()
        return li
