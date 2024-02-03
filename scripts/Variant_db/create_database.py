#!/usr/bin/python3
import sys
import sqlite3

def create_database():
    if len(sys.argv) != 2:
        print("Usage: {} <database_name>".format(sys.argv[0]))
        sys.exit(1)

    db_name = sys.argv[1]
    conn = sqlite3.connect(db_name)  # Removed the $ sign
    cursor = conn.cursor()

    # Create the Variants table without the GT column
    cursor.execute('''
        CREATE TABLE IF NOT EXISTS Variants (
            CHROM TEXT,
            POS INTEGER,
            REF TEXT,
            ALT TEXT,
            SAMPLE_NAME TEXT,
            ALT_NUM INTEGER,
            PRIMARY KEY (CHROM, POS, REF, ALT)
        );
    ''')

    conn.commit()
    conn.close()

if __name__ == "__main__":
    create_database()
    print("Database '{}' created with 'Variants' table".format(sys.argv[1]))

