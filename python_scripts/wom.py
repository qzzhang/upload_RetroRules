
import sqlite3
from sqlite3 import Error
import csv
import pandas as pd
import datetime
import pprint


def create_connection(db_file):
    """ create a database connection to the SQLite database
        specified by the db_file
    :param db_file: database file
    :return: Connection object or None
    """
    try:
        conn = sqlite3.connect(db_file)
        return conn
    except Error as e:
        print(e)
 
    return None


def column_attribute_query():
    """
    column_attribute_query: construct a query that retrieves data for defining
    the matrix's column attributes from four tables matchmaker_observation
    /organism/environment/project
    :return: an sqlite3 query string
    """
    qry = """
    select distinct 'eop_E'||ob.environment_id||'-O'||ob.organism_id||'-P'||ob.project_id as env_org_proj_id,
    en.env_name as environment_name, org.common_name as organism_name, prj.project_name as project_name,
    org.NCBI_taxid, prj.contributor, prj.project_description
    from matchmaker_observation ob,
    matchmaker_environment en,
    matchmaker_project prj,
    matchmaker_organism org
    where ob.organism_id=org.id and ob.project_id=prj.id and ob.environment_id=en.id
    """
    return qry


def row_attribute_query():
    """
    row_attribute_query: construct a query that retrieves data for defining
    the matrix's row attributes from one table: matchmaker_compound
    :return: an sqlite3 query string
    """
    qry = """
    SELECT compound_name as cpd_name, formula, id as cpd_id, neutralmass
    from matchmaker_compound
    """
    # where formula!='NA' and compound_name not like 'Unk_%'
    return qry


def matrix_query():
    """
    matrix_query: construct a query that retrieves data for defining
    the matrix's values from four tables matchmaker_observation
    /organism/environment/project
    :return: an sqlite3 query string
    """
    qry = """
        select group_concat(compound_id) as cpd_ids,
        group_concat(action) as actions, group_concat(confidence) as confidences,
        'eop_E'||environment_id||'-O'||organism_id||'-P'||project_id as env_org_proj_id
        from matchmaker_observation
        group by env_org_proj_id
    """
    return qry


def post_query_process(data_rows, row_count=0):
    """
    post_query_process: Further massage the data to meet with
    downstream input required formats
    :param data_rows : SQL query result in a format of list of tupples
    :param row_count : The number of rows out of data_rows to be processed
    :return: a dict object and a list object
    """
    in_data = []
    dict_data = {}
    list_data = []
    if row_count > 0:
        in_data = data_rows[:row_count]
    else:
        in_data = data_rows

    # pre-fill the first column using 'compound_id' as caption and
    # compound_id's from 1 through row_count.
    dict_data['compound_id'] = [i + 1 for i in range(row_count)]

    cpd_ids = []
    actions = []
    confidences = []
    for row in in_data:
        lst_row = list(row)
        out_data_col = [None] * row_count
        cpd_ids = lst_row[0].split(',')
        actions = lst_row[1].split(',')
        confidences = lst_row[2].split(',')

        for i in range(len(cpd_ids)):
            int_cpd_id = int(cpd_ids[i])
            if actions[i] == 'N':
                # Sample metabolite not detected or no significant change
                out_data_col[int_cpd_id - 1] = ''
            elif actions[i] == 'D':  # Control metabolite detected
                out_data_col[int_cpd_id - 1] = ''
            elif actions[i] == 'I':  # Intake
                out_data_col[int_cpd_id - 1] = str(-float(confidences[i]))
            elif actions[i] == 'E':  # Excrete
                out_data_col[int_cpd_id - 1] = str(float(confidences[i]))
        k_val = lst_row[3]
        if ',' in lst_row[3]:
            k_val = '"{}"'.format(k_val)
        dict_data[lst_row[3]] = out_data_col
        list_data.append(out_data_col)
    return (dict_data, list_data)


def process_col_results(data_rows, row_count=0):
    """
    process_col_results: Further massage the data to meet with
    downstream input required formats
    :param data_rows : SQL query result in a format of list of tupples
    :param row_count : The number of rows out of data_rows to be processed
    :return: a dict object and a list object
    """
    in_data = []
    out_data = []
    if row_count > 0:
        in_data = data_rows[:row_count]
    else:
        in_data = data_rows

    for row in in_data:
        lst_row = list(row)
        if ',' in lst_row[0]:
            lst_row[0] = '"{}"'.format(lst_row[0])
        out_data.append(tuple(lst_row))

    return out_data


def execute_query(conn, qry):
    """
    execute_query: executes the given query qry against db connection con
    :param conn: the Connection object
    :param qry: SQL query smarts_string
    :return: a list of rows (tuples) if no error, otherwise None
    """
    ret_data = None
    # if sqlite3.complete_statement(qry):
    try:
        cur = conn.cursor()
        qry = qry.strip()
        cur.execute(qry)
        ret_data = cur.fetchall()
    except Error as e:
        print("An error occurred:", e.args[0])

    return ret_data


def build_matrix(conn, row_count):
    """
    build_matrix: construct a matrix that matches data from querying
    table matchmaker_observation with conditions in the given
    matrix's column attributes cols
    :return: a 2-D array
    """
    qry = matrix_query()
    return post_query_process(execute_query(conn, qry), row_count)


def row_attribute_process(in_data, cpd_inchi_dict):
    """
    row_attribute_process: insert the inchiKey column to the input in_data
    :param in_data : SQL query result in a format of list of tupples
    :param cpd_inchi_dict : The dict for lookup with a cpd name for inchiKey
    """
    out_data = []
    for row in in_data:
        lst_row = list(row)
        lst_row.append(cpd_inchi_dict.get(lst_row[0], ''))
        out_data.append(lst_row)
    return out_data


def csv_dict_reader(file_obj):
    """
    Read a CSV file using csv.DictReader
    We open a file and pass the file object to our function as we did before. The function passes
    the file object to our DictReader class. We tell the DictReader that the delimiter is a comma.
    This isn't actually required as the code will still work without that keyword argument. However,
    it's a good idea to be explicit so you know what's going on here. Next we loop over the reader
    object and discover that each line in the reader object is a dictionary. This makes printing out
    specific pieces of the line very easy.
    """
    reader = csv.DictReader(file_obj, delimiter=',')
    for line in reader:
        print(line["first_name"]),
        print(line["last_name"])


def read_tsv_into_dict(tsvfile_path):
    """
    read_tsv_into_dict: read info from given tsv file path
    :param tsvfile_path: filename with path to read from
    :return: a dictionary of compound to inchiKey
    """
    cpd_inchis = {}
    with open(tsvfile_path, "r") as file_obj:
        reader = csv.reader(file_obj)
        for row in reader:
            # print(row)  # each row is an array of a single tab delimited string
            row_arr = row[0].split("\t")
            cpd_inchis[row_arr[0]] = row_arr[2]
    # print(cpd_inchis)  # .get("valine"))
    return cpd_inchis


def csv_write(data, fpath, fheader):
    """
    csv_write: Write data to a CSV file path
    :param data: a list of lists that is returned by a database query
    :param fpath: filename with path to write to
    :param fheader: data header in the TSV file
    e.g., ["cpd_name", "formula", "cpd_id", "neutralmass", "inchiKey"]
    """
    with open(fpath, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter='\t')  # create a csv.writer, tab delimited
        # write the header
        writer.writerow(fheader)
        # write the data rows
        writer.writerows(row for row in data)


def csv_dict_writer(path, fieldnames, data, delm):
    """
    csv_dict_writer: Writes a CSV file using DictWriter
    """
    with open(path, "wb") as out_file:
        # create a csv.DictWriter
        writer = csv.DictWriter(out_file, delimiter=delm,
                                fieldnames=fieldnames)
        writer.writeheader()

        # Here we could use writer.writerrows() instead of the loop
        # i.e.,
        # enc_data = dict((k, v.encode('utf-8') if isinstance(v, unicode) else v)
        #                  for k, v in data.iteritems())
        # writer.writerows(enc_data)
        for row in data:
            writer.writerow(row)


def main():
    # pp = pprint.PrettyPrinter(indent=4)
    database = "/Users/qzhang/qzwk_dir/wom/wom.sqlite3"
    print("1. create a database connection...")
    conn = create_connection(database)

    with conn:
        now = datetime.datetime.now()
        date_str = now.strftime("%Y-%m-%d")

        print("2.1 Column query result...")
        col_qry_result = execute_query(conn, column_attribute_query())
        eop_header = ["env_org_proj_id", "environment_name", "organism_name",
                      "project_name", "NCBI_taxid", "contributor",
                      "project_description"]
        eopfile_nm = "../TSVs/wom_eop_{}{}".format(date_str, ".tsv")
        print("2.2 Write environment_organism_project info to output file {}"
              .format(eopfile_nm))
        if col_qry_result:
            # pp.pprint(col_qry_result)
            csv_write(col_qry_result, eopfile_nm, eop_header)

        print("3.1. read inchikey file...")
        csvfile_path = "/Users/qzhang/qzwk_dir/wom/genericLoading/" + \
            "formula_withoutNA_withInchiKeys.tsv"

        cpd_inchi_dict = read_tsv_into_dict(csvfile_path)

        print("3.2 Row query result...")
        row_qry_result = execute_query(conn, row_attribute_query())

        print("3.3 Inserting inchikey to compounds")
        row_result = row_attribute_process(row_qry_result, cpd_inchi_dict)
        cpd_header = ["cpd_name", "formula", "cpd_id",
                      "mass", "inchikey"]
        cpdfile_nm = "../TSVs/wom_cpd_{}{}".format(date_str, ".tsv")

        print("3.4 Write compounds to output file {}".format(cpdfile_nm))
        if row_result:
            csv_write(row_result, cpdfile_nm, cpd_header)

        row_count = len(row_result)
        matrix_dict, matrix_list = build_matrix(conn, row_count)
        matrixfile_nm = "../TSVs/wom_matrix_df_{}{}".format(date_str, ".tsv")

        print("4. Write matrix to output file {}".format(matrixfile_nm))
        df = pd.DataFrame(matrix_dict, columns=list(matrix_dict.keys()))
        df.to_csv(matrixfile_nm, sep='\t', index=False, decimal='.')


if __name__ == '__main__':
    """
    Note: After generate the wom_cpd_*.tsv/wom_eop_*.tsv/wom_matrix_df_*.tsv
          files, there is still another step to do before uploading the tsv
          files to KBase.Namely, the first columns of wom_cpd_*.tsv and
          wom_matrix_df_*.tsv have tomatch with each other value-wise,
          even though their column caption names can be different
          (e.g., compound_name or compound_id column for both or
          for either file).
    """
    main()
    