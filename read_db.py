import sqlite3

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

def read_residue_db(selection, key, value):
    connection = sqlite3.connect("/Users/mscheurer/Projects/pycontact/pycontact/aa.db")
    connection.row_factory = dict_factory
    cursor = connection.cursor()
    cursor.execute("SELECT ({sel}) FROM residues WHERE {key}=\'{value}\'".format(sel=selection, key=key, value=value))
    results = cursor.fetchall()
    connection.close()
    return results