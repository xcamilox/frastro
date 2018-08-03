from astropy.io.votable.tree import VOTableFile, Resource, Table, Field
from astropy.io.votable import parse_single_table
from astropy.io.votable import from_table, writeto

class VOTableUtil():
    @staticmethod
    def createTableFromObject(data,path="",names=[],dtypes=[],sizes=[]):

        path_tmp = "/Users/cjimenez/Documents/PHD/data/tmp/"
        # Create a new VOTable file...
        votable = VOTableFile()

        # ...with one resource...
        resource = Resource()
        votable.resources.append(resource)

        # ... with one table
        table = Table(votable)
        resource.tables.append(table)

        # Define some fields

        fields = []
        for idx, val in enumerate(names):
            fields.append(Field(votable, name=val, datatype=dtypes[idx]))

        table.fields.extend(fields)

        # Now, use those field definitions to create the numpy record arrays, with
        # the given number of rows
        table.create_arrays(len(data))

        # Now table.array can be filled with data
        for idx, val in enumerate(data):
            table.array[idx] = val


        # Now write the whole thing to a file.
        # Note, we have to use the top-level votable file object
        votable.to_xml(path_tmp+path)

    @staticmethod
    def saveFromTable(table,ouputfile):
        votable = from_table(table)
        writeto(votable, ouputfile)