# https://docs.python.org/2/library/sqlite3.html#sqlite3.Connection.row_factory

import sqlite3

def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d

connection = sqlite3.connect("sample.db")
connection.row_factory = dict_factory
cursor = connection.cursor()
cursor.execute("SELECT name FROM residues WHERE hbondtype='none'")
results = cursor.fetchall()
print(results)
connection.close()
