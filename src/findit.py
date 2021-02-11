def search(dic, searchFor):
    for key, value in dic.items():
        if isinstance(value, str):
            if searchFor in value: print(key,value)
                
    return None

from distutils.sysconfig import *
dic = get_config_vars()

search(dic,'lib')
print(dic['CFLAGS'])
      
