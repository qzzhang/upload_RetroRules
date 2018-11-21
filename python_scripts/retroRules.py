
import sqlite3
from sqlite3 import Error
import itertools
import csv
import re


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


def build_query(diam=10):
    """
    build_query: construct a query across tables rules, rule_products, reaction, smarts,
    reaction_substrates, reaction_products, chemical_species and ec_numbers
    :param diam: reaction diameter
    :return: an sqlite3 query string
    """
    rl_info = """
                select rl.reaction_id,rl.substrate_id,rl.rule_substrate_cpd,rl.diameter,rl.direction,rl.isStereo,rl.score,rl.SMARTS,
                product_per_rxn_sub_dia_isStereo.rule_prod_ids,product_per_rxn_sub_dia_isStereo.rule_prod_stoichios,
                product_per_rxn_sub_dia_isStereo.total_stoichios
                from
                (
                    select r.reaction_id,r.substrate_id,r.score,r.isStereo,s.smarts_string as SMARTS,r.diameter,r.direction,
                    cs_info.repo_cpd_id as rule_substrate_cpd
                    from rules r, smarts s,
                    (
                        select cs.id, ifnull(cs.seed, ifnull(cs.bigg, ifnull(cs.kegg, ifnull(cs.metacyc, cs.mnxm)))) as repo_cpd_id
                        from chemical_species cs
                    ) as cs_info
                    where s.id=r.smarts_id and cs_info.id=r.substrate_id
                ) as rl,
                (
                    select reaction_id,substrate_id,diameter,isStereo,
                    group_concat(repo_prod_id) as rule_prod_ids,
                    group_concat(prod_stoichio) as rule_prod_stoichios,
                    sum(prod_stoichio) as total_stoichios
                    from
                    (
                        select rp.reaction_id, rp.substrate_id, rp.diameter,rp.isStereo,
                        cs_info1.repo_cpd_id as repo_prod_id, rp.stochiometry as prod_stoichio
                        from rule_products rp,
                        (
                            select cs.id, ifnull(cs.seed, ifnull(cs.bigg, ifnull(cs.kegg, ifnull(cs.metacyc, cs.mnxm)))) as repo_cpd_id
                            from chemical_species cs
                        ) as cs_info1
                        where rp.product_id=cs_info1.id
                    )
                    group by reaction_id,substrate_id,diameter,isStereo
                ) as product_per_rxn_sub_dia_isStereo
                where (rl.reaction_id=product_per_rxn_sub_dia_isStereo.reaction_id and
                rl.substrate_id=product_per_rxn_sub_dia_isStereo.substrate_id and
                rl.diameter=product_per_rxn_sub_dia_isStereo.diameter and
                rl.isStereo=product_per_rxn_sub_dia_isStereo.isStereo)
    """
    rl_info1 = """(
        select rl_info.reaction_id,rl_info.substrate_id as rule_substrate_id,rl_info.rule_substrate_cpd,rl_info.direction,
        rl_info.diameter,rl_info.isStereo,rl_info.score,rl_info.SMARTS,rl_info.rule_prod_ids,
        rl_info.rule_prod_stoichios,rl_info.total_stoichios
        from (""" + rl_info + ") as rl_info where rl_info.diameter=" + str(diam) + ") as rl_info1"

    rxn_info = """(
        select rxn_substrates_tab.id,rxn_substrates_tab.repo_rxn_id,rxn_substrates_tab.ec_numbers,
        rxn_substrates_tab.substrate_ids as rxn_substrate_ids,rxn_substrates_tab.substrate_inchi_keys as rxn_substrate_inchis,
        rxn_products_tab.product_ids as rxn_product_ids,rxn_products_tab.product_inchi_keys as rxn_product_inchis
        from
        (
            select rxn.id, rxn.repo_rxn_id,group_concat(distinct rxn.ec_number) as ec_numbers,
            group_concat(distinct rxn_substrates.substrate_cpd) as substrate_ids,
            group_concat(distinct rxn_substrates.inchi_key) as substrate_inchi_keys
            from
            (
                select rxn2.id, rxn2.repo_rxn_id,er.ec_number
                from
                (
                    select rxn1.id, ifnull(rxn1.seed, ifnull(rxn1.bigg, ifnull(rxn1.kegg, ifnull(rxn1.metacyc, rxn1.mnxr)))) as repo_rxn_id
                    from reactions rxn1
                ) as rxn2
                left join ec_reactions er
                on er.reaction_id=rxn2.id
            ) as rxn
            left join
            (
                select rs.reaction_id, cs_info.cpd_id as substrate_cpd, rs.chemical_id, cs_info.inchi_key
                from
                (
                    select distinct cs.id as chem_sp_id, cs.inchi_key,ifnull(cs.seed, ifnull(cs.bigg, ifnull(cs.kegg, ifnull(cs.metacyc, cs.mnxm)))) as cpd_id
                    from chemical_species cs
                ) as cs_info, reaction_substrates rs
                where rs.chemical_id=cs_info.chem_sp_id
            ) as rxn_substrates
            on rxn.id=rxn_substrates.reaction_id
            group by rxn.id
        ) as rxn_substrates_tab,
        (
            select rxn.id, rxn.repo_rxn_id,group_concat(distinct rxn.ec_number) as ec_numbers,
            group_concat(distinct rxn_products.product_cpd) as product_ids,
            group_concat(distinct rxn_products.inchi_key) as product_inchi_keys
            from
            (
                select rxn2.id, rxn2.repo_rxn_id,er.ec_number
                from
                (
                    select rxn1.id, ifnull(rxn1.seed, ifnull(rxn1.bigg, ifnull(rxn1.kegg, ifnull(rxn1.metacyc, rxn1.mnxr)))) as repo_rxn_id
                    from reactions rxn1
                ) as rxn2
                left join ec_reactions er
                on er.reaction_id=rxn2.id
            ) as rxn
            left join
            (
                select rp.reaction_id, cs_info.cpd_id as product_cpd,rp.chemical_id,cs_info.inchi_key
                from
                (
                    select distinct cs.id as chem_sp_id, cs.inchi_key,ifnull(cs.seed, ifnull(cs.bigg, ifnull(cs.kegg, ifnull(cs.metacyc, cs.mnxm)))) as cpd_id
                    from chemical_species cs
                ) as cs_info, reaction_products rp
                where rp.chemical_id=cs_info.chem_sp_id
            ) as rxn_products
            on rxn.id=rxn_products.reaction_id
            group by rxn.id
        ) as rxn_products_tab
        where rxn_substrates_tab.id=rxn_products_tab.id
    ) as rxn_info
    """

    qry = rl_info1 + " left join " + rxn_info + " on rxn_info.id=rl_info1.reaction_id "

    qry = """
    select rl_info1.reaction_id,rxn_info.repo_rxn_id,rxn_info.ec_numbers,
        rl_info1.rule_substrate_id,rl_info1.rule_substrate_cpd,
        rl_info1.direction,rl_info1.diameter,rl_info1.isStereo,rl_info1.score,
        rl_info1.SMARTS, 'Any' as Reactants, rl_info1.rule_prod_ids,
        rl_info1.rule_prod_stoichios,rxn_info.rxn_substrate_ids,rxn_info.rxn_substrate_inchis,
        rxn_info.rxn_product_ids,rxn_info.rxn_product_inchis,rl_info1.total_stoichios,
        (select '<'||rl_info1.reaction_id||'|'||rxn_info.repo_rxn_id||'|'||
            rl_info1.rule_substrate_cpd||'|'||rl_info1.diameter||'|'||
            (SELECT
            CASE
                WHEN rl_info1.direction=-1 THEN 'reverse'
                ELSE 'forward'
            END) ||
            (SELECT
            CASE
                WHEN rl_info1.isStereo=1 THEN '|'||'isStereo'||'>'
                ELSE ''||'>'
            END) as Name)
    from """ + qry

    return qry


def build_query_seed_cpds(diam=10):
    """
    build_query_seed_cpds: construct a query across tables rules, rule_products, reaction, smarts,
    reaction_substrates, reaction_products, chemical_species and ec_numbers ONLY for rules that
    have seed reactants/products
    :param diam: reaction diameter
    :return: an sqlite3 query string
    """
    rl_info = """
                select rl.reaction_id,rl.substrate_id,rl.rule_substrate_cpd,rl.diameter,rl.direction,rl.isStereo,rl.score,rl.SMARTS,
                product_per_rxn_sub_dia_isStereo.rule_prod_ids,product_per_rxn_sub_dia_isStereo.rule_prod_stoichios,
                product_per_rxn_sub_dia_isStereo.total_stoichios
                from
                (
                    select r.reaction_id,r.substrate_id,r.score,r.isStereo,s.smarts_string as SMARTS,r.diameter,r.direction,
                    cs_info.repo_cpd_id as rule_substrate_cpd
                    from rules r, smarts s,
                    (
                        select cs.id, cs.seed as repo_cpd_id
                        from chemical_species cs
                        where cs.seed not null
                    ) as cs_info
                    where s.id=r.smarts_id and cs_info.id=r.substrate_id
                ) as rl,
                (
                    select reaction_id,substrate_id,diameter,isStereo,
                    group_concat(repo_prod_id) as rule_prod_ids,
                    group_concat(prod_stoichio) as rule_prod_stoichios,
                    sum(prod_stoichio) as total_stoichios
                    from
                    (
                        select rp.reaction_id, rp.substrate_id, rp.diameter,rp.isStereo,
                        cs_info1.repo_cpd_id as repo_prod_id, rp.stochiometry as prod_stoichio
                        from rule_products rp,
                        (
                            select cs.id, cs.seed as repo_cpd_id
                            from chemical_species cs
                            where cs.seed not null
                        ) as cs_info1
                        where rp.product_id=cs_info1.id
                    )
                    group by reaction_id,substrate_id,diameter,isStereo
                ) as product_per_rxn_sub_dia_isStereo
                where (rl.reaction_id=product_per_rxn_sub_dia_isStereo.reaction_id and
                rl.substrate_id=product_per_rxn_sub_dia_isStereo.substrate_id and
                rl.diameter=product_per_rxn_sub_dia_isStereo.diameter and
                rl.isStereo=product_per_rxn_sub_dia_isStereo.isStereo)
    """
    rl_info1 = """(
        select rl_info.reaction_id,rl_info.substrate_id as rule_substrate_id,rl_info.rule_substrate_cpd,rl_info.direction,
        rl_info.diameter,rl_info.isStereo,rl_info.score,rl_info.SMARTS,rl_info.rule_prod_ids,
        rl_info.rule_prod_stoichios,rl_info.total_stoichios
        from (""" + rl_info + ") as rl_info where rl_info.diameter=" + str(diam) + ") as rl_info1"

    rxn_info = """(
        select rxn_substrates_tab.id,rxn_substrates_tab.repo_rxn_id,rxn_substrates_tab.ec_numbers,
        rxn_substrates_tab.substrate_ids as rxn_substrate_ids,rxn_substrates_tab.substrate_inchi_keys as rxn_substrate_inchis,
        rxn_products_tab.product_ids as rxn_product_ids,rxn_products_tab.product_inchi_keys as rxn_product_inchis
        from
        (
            select rxn.id, rxn.repo_rxn_id,group_concat(distinct rxn.ec_number) as ec_numbers,
            group_concat(distinct rxn_substrates.substrate_cpd) as substrate_ids,
            group_concat(distinct rxn_substrates.inchi_key) as substrate_inchi_keys
            from
            (
                select rxn2.id, rxn2.repo_rxn_id,er.ec_number
                from
                (
                    select rxn1.id,  ifnull(rxn1.seed, ifnull(rxn1.bigg, ifnull(rxn1.kegg, ifnull(rxn1.metacyc, rxn1.mnxr)))) as repo_rxn_id
                    from reactions rxn1
                ) as rxn2
                left join ec_reactions er
                on er.reaction_id=rxn2.id
            ) as rxn
            left join
            (
                select rs.reaction_id, cs_info.cpd_id as substrate_cpd, rs.chemical_id, cs_info.inchi_key
                from
                (
                    select distinct cs.id as chem_sp_id, cs.inchi_key, cs.seed as cpd_id
                    from chemical_species cs
                    where cs.seed not null
                ) as cs_info, reaction_substrates rs
                where rs.chemical_id=cs_info.chem_sp_id
            ) as rxn_substrates
            on rxn.id=rxn_substrates.reaction_id
            group by rxn.id
        ) as rxn_substrates_tab,
        (
            select rxn.id, rxn.repo_rxn_id,group_concat(distinct rxn.ec_number) as ec_numbers,
            group_concat(distinct rxn_products.product_cpd) as product_ids,
            group_concat(distinct rxn_products.inchi_key) as product_inchi_keys
            from
            (
                select rxn2.id, rxn2.repo_rxn_id,er.ec_number
                from
                (
                    select rxn1.id, ifnull(rxn1.seed, ifnull(rxn1.bigg, ifnull(rxn1.kegg, ifnull(rxn1.metacyc, rxn1.mnxr)))) as repo_rxn_id
                    from reactions rxn1
                ) as rxn2
                left join ec_reactions er
                on er.reaction_id=rxn2.id
            ) as rxn
            left join
            (
                select rp.reaction_id, cs_info.cpd_id as product_cpd, rp.chemical_id, cs_info.inchi_key
                from
                (
                    select distinct cs.id as chem_sp_id, cs.inchi_key, cs.seed as cpd_id
                    from chemical_species cs
                    where cs.seed not null
                ) as cs_info, reaction_products rp
                where rp.chemical_id=cs_info.chem_sp_id
            ) as rxn_products
            on rxn.id=rxn_products.reaction_id
            group by rxn.id
        ) as rxn_products_tab
        where rxn_substrates_tab.id=rxn_products_tab.id
    ) as rxn_info
    """

    qry = rl_info1 + " left join " + rxn_info + " on rxn_info.id=rl_info1.reaction_id "

    qry = """
    select rl_info1.reaction_id,rxn_info.repo_rxn_id,rxn_info.ec_numbers,
        rl_info1.rule_substrate_id,rl_info1.rule_substrate_cpd,
        rl_info1.direction,rl_info1.diameter,rl_info1.isStereo,rl_info1.score,
        rl_info1.SMARTS, 'Any' as Reactants, rl_info1.rule_prod_ids,
        rl_info1.rule_prod_stoichios,rxn_info.rxn_substrate_ids,rxn_info.rxn_substrate_inchis,
        rxn_info.rxn_product_ids,rxn_info.rxn_product_inchis,rl_info1.total_stoichios,
        (select '<'||rl_info1.reaction_id||'|'||rxn_info.repo_rxn_id||'|'||
            rl_info1.rule_substrate_cpd||'|'||rl_info1.diameter||'|'||
            (SELECT
            CASE
                WHEN rl_info1.direction=-1 THEN 'reverse'
                ELSE 'forward'
            END) ||
            (SELECT
            CASE
                WHEN rl_info1.isStereo=1 THEN '|'||'isStereo'||'>'
                ELSE ''||'>'
            END) as Name)
    from """ + qry

    return qry


def execute_query(conn, qry):
    """
    execute_query: executes the given query qry against db connection conn
    :param conn: the Connection object
    :param qry: SQL query smarts_string
    :return: a list of rows (tuples) if no error, otherwise None
    """
    ret_data = None
    # if sqlite3.complete_statement(qry):
    # print(qry)
    try:
        cur = conn.cursor()
        qry = qry.strip()
        cur.execute(qry)
        rows = cur.fetchall()
        ret_data = rows
    except Error as e:
        print("An error occurred:", e.args[0])

    return ret_data


def repeat_any(N):
    """
    repeat_Any: return a string with 'Any' repeated N times
    """
    rep_str = 'Any'  # Start with one 'Any'
    for _ in itertools.repeat(None, N-1):
        rep_str += ';Any'
    return rep_str


def post_query_process(data_rows, row_count=0):
    """
    postQueryProcess: Further massage the data to meet with
    downstream input required formats
    :param data_rows : SQL query result in a format of list of tupples
    :param row_count : The number of rows out of in_data to be processed
    """
    in_data = []
    out_data = []
    if row_count > 0:
        in_data = data_rows[:row_count]
    else:
        in_data = data_rows

    pattern = r'>>\((.*)\)$'
    for row in in_data:
        lst_row = list(row)
        any_num = lst_row[17]  # 'total_stoichios'
        smt = lst_row[9]  # 'SMARTS'
        if any_num > 1:
            lst_row[9] = re.sub(pattern, r'>>\1', smt)
        lst_row[17] = (repeat_any(any_num))
        out_data.append(lst_row)

    return out_data


def generate_rule_per_row_table(conn, row_count=0, diam=10):
    """
    generate_rule_per_row_table: Query the tables rules, rule_products, reactions,
    reaction_substrates, reaction_products, smarts, chemical_species and ec_numbers
    :param conn: the Connection object
    :param row_count: number of rows to output, if 0 return all
    :param diam: reaction diameter
    :return: a list of rows (tuples) if no error, otherwise None
    """
    qry = build_query(diam)
    # qry_seed = (qry + " where rxn_info.repo_rxn_id='rxn14222'" +
    #           " and rl_info1.rule_substrate_cpd='cpd17740'")
    qry_seed = qry + " where rxn_info.repo_rxn_id like 'rxn1%'"

    qry_result = execute_query(conn, qry_seed)

    return post_query_process(qry_result, row_count)


def generate_rule_per_row_table_seed_cpds(conn, row_count=0, diam=10):
    """
    generate_rule_per_row_table: Query the tables rules, rule_products, reactions,
    reaction_substrates, reaction_products, smarts, chemical_species and ec_numbers
    :param conn: the Connection object
    :param row_count: number of rows to output, if 0 return all
    :param diam: reaction diameter
    :return: a list of rows (tuples) if no error, otherwise None
    """
    qry = build_query_seed_cpds(diam)
    qry_result = execute_query(conn, qry)

    return post_query_process(qry_result, row_count)


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
        print(line["first_name"])
        print(line["last_name"])


def csv_write(data, fpath):
    """
    Write data to a CSV file path
    We create a csv_writer function that accepts two arguments: data and path.
    The data is a list of lists that is returned by a database query.
    :param fpath: filename with path to write to
    """
    with open(fpath, "w") as csv_file:
        writer = csv.writer(csv_file, delimiter='\t')  # create a csv.writer, tab delimited
        # write the header
        writer.writerow(
            ["reaction_id", "repo_rxn_id", "ec_numbers", "rule_substrate_id", "rule_substrate_cpd",
             "direction", "diameter", "isStereo", "score", "SMARTS", "Reactants", "rule_prod_ids",
             "rule_prod_stoichios", "rxn_substrate_ids", "rxn_substrate_inchis", "rxn_product_ids",
             "rxn_product_inchis", "Products", "Name"])
        # write the data rows
        writer.writerows(row for row in data)


def csv_dict_writer(path, fieldnames, data, delm):
    """
    Writes a CSV file using DictWriter
    """
    with open(path, "wb") as out_file:
        # create a csv.DictWriter
        writer = csv.DictWriter(out_file, delimiter=delm,
                                fieldnames=fieldnames)
        writer.writeheader()

        # Here we could use writer.writerrows() instead of the loop
        # i.e.,
        # enc_data = dict((k, v.encode('utf-8') if isinstance(v, unicode) else v) for k, v in data.iteritems())
        # writer.writerows(enc_data)
        for row in data:
            writer.writerow(row)


def main():
    """
    main: The main
    """
    database = "/Users/qzhang/qzwk_dir/upload_RetroRules/retrorules_dump/mvc.db"

    print("0. create a database connection...")
    conn = create_connection(database)
    row_cnt = 0   # 200
    diam = 16
    str_row_cnt = str(row_cnt) if row_cnt > 0 else 'all'
    outfile_nm = "../TSVs/retro_rules_dia{}_{}".format(str(diam), str_row_cnt) + ".tsv"
    outfile_nm_seed_cpds = "../TSVs/seed_cpds_rules_dia{}_{}".format(str(diam), str_row_cnt) + ".tsv"

    rule_per_row_results = None
    rule_per_row_results_seed_cpds = None
    with conn:
        print("1. Query tables to create the result data...")
        # rule_per_row_results = generate_rule_per_row_table(conn, row_cnt, diam)
        rule_per_row_results_seed_cpds = generate_rule_per_row_table_seed_cpds(conn, row_cnt, diam)
            
        print("2. Write to output file {}".format(outfile_nm))
        if rule_per_row_results:
            csv_write(rule_per_row_results, outfile_nm)

        print("2. Write to output file {}".format(outfile_nm_seed_cpds))
        if rule_per_row_results_seed_cpds:
            csv_write(rule_per_row_results_seed_cpds, outfile_nm_seed_cpds)


if __name__ == '__main__':
    main()