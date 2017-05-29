import sqlite3
import os


def dict_factory(cursor, row):
    d = {}
    for idx, col in enumerate(cursor.description):
        d[col[0]] = row[idx]
    return d


def read_residue_db(selection, key, value):
    connection = sqlite3.connect(os.path.dirname(os.path.abspath(__file__)) + "/aa.db")
    connection.row_factory = dict_factory
    cursor = connection.cursor()
    cursor.execute("SELECT ({sel}) FROM residues WHERE {key}=\'{value}\'".format(sel=selection, key=key, value=value))
    results = cursor.fetchall()
    connection.close()
    return results


def read_residue_db_all():
    connection = sqlite3.connect(os.path.dirname(os.path.abspath(__file__)) + "/aa.db")
    connection.row_factory = dict_factory
    cursor = connection.cursor()
    cursor.execute("SELECT * FROM residues")
    results = cursor.fetchall()
    connection.close()
    return results
