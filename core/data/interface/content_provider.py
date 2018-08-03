import abc

"""
Content Provider
"""
class ContentProvider(metaclass=abc.ABCMeta):

    @property
    @abc.abstractmethod
    def contentUrl(self):
        return str

    @property
    @abc.abstractmethod
    def getType(self):
        pass

    @abc.abstractmethod
    def onCreate(self):
        pass

    @abc.abstractmethod
    def query(self, **kwargs):
        pass

    @abc.abstractmethod
    def insert(self, **kwargs):
        pass

    @abc.abstractmethod
    def update(self, **kwargs):
        pass

    @abc.abstractmethod
    def delete(self, **kwargs):
        pass


