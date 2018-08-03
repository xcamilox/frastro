from bs4 import BeautifulSoup
class HtmlParser():
    @staticmethod
    def getImageList(html_str):
        return HtmlParser.findAllElements(html_str, "img")

    @staticmethod
    def findAllElements(html_str,element_str):
        parser = BeautifulSoup(html_str)
        list = parser.find_all(element_str)
        return list

    @staticmethod
    def getLink(html_str):
        return HtmlParser.findAllElements(html_str,"a")

    @staticmethod
    def getSource(html_str):
        return HtmlParser.findAllElements(html_str, "src")

    @staticmethod
    def getById(html_str,id_str):
        parser = BeautifulSoup(html_str,'html.parser')
        list = parser.find(id=id_str)
        return list

