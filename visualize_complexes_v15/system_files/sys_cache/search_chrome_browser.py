import os


def get_bowser_location():

	main_directory = 'C:\\'
	browser_exe = 'chrome.exe'

	user = os.getlogin()
	lis_dir = ['C:/Program Files (x86)/Google/Chrome/Application',
			   'C:/Program Files/Google/Chrome/Application',
			   'C:/Users/' + user + '/AppData/Local/Google/Chrome/Application']

	def search_usual():
		for dir in lis_dir:
			try:
				for file in os.listdir(dir):
					if file.endswith(browser_exe):
						return dir + '/' + file
			except FileNotFoundError:
				return None

	file_loc = search_usual()

	return file_loc
