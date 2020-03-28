import settings, os, sys

def test_settings():
	for k in settings.LOC.keys():
		print(k)
		assert os.path.isdir(settings.LOC[k])
	return

if __name__=='__main__':
	test_settings()
