# Dirty hack to make nglview 3.1.4 work in a Python distribution without pkg_resources
class Distribution(object):
    def __init__(self, version=None):                
        self.version = version
        
def get_distribution(dist):
    dist = Distribution(version="3.1.4")
    return dist
    