import pdfkit
class ReportUtil():
    html = """
           <html>
             <head>
               <meta name="pdfkit-page-size" content="Legal"/>
               <meta name="pdfkit-orientation" content="Landscape"/>
             </head>
             <body>
         """

    def __init__(self,url_report=""):
        self.output = url_report

    def saveReport(self,output_path,content=""):

        if content == "":
            html = self.html
        else:
            html = content

        html += """
                 </body>
             </html>
           """
        f = open(output_path, "w+")
        f.write(html)
        f.close()


    def addImage(self,img_url,width=100):
        return "<img src='" + img_url + "' width='"+str(width)+"'>"

    def addText(self,text):
        return "<p>"+text+"</p>"

    def addTitle(self,title):
        return "<h1>"+title+"</h1>"

    def addSubTitle(self, subTitle):
        return "<h3>" + subTitle + "</h3>"

    def wrapContent(self,content):
        return "<div>"+content+"</div>"

    def createPDF(self,file_path,str_content=""):
        if str_content == "":
            html_str = self.html
        else:
            html_str = str_content
        pdfkit.from_string(html_str, file_path)

    def createTable(self,header=[]):
        table = '<table style="width:100%"><thead><tr>'
        header_title=""
        for title in header:
            header_title+="<th>"+str(title)+"</th>"
        table += header_title + "</tr></thead><tbody>"
        return table

    def addToTable(self,rows=[],img=""):
        row="<tr>"
        items=""
        for content in rows:
            items+="<td>"+str(content)+"</td>"

        if img!="":
            items += '<td><img src='+ str(img) +' height="50%"></td>'

        row+=items+"</tr>"
        return row

    def saveTable(self,table):
        table+="</tbody>"
        return table


