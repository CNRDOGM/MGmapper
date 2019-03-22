#!/usr/bin/env python2.7
# Created: Thu Mar  6 19:50:12 2014
# Last changed: Time-stamp: <2014-10-13 12:52:55 tnp>
# Latest patch: Line 49 - unicode characters <2015-10-13 11:16 koala>
# thomas@cbs.dtu.dk, http://www.cbs.dtu.dk/thomas
import os, sys, commands, re
import argparse
#sys.path.insert(0, '/usr/local/python/lib')
sys.path.insert(0, os.environ['MGmap_HOME'] + '/modules')
import xlsxwriter

def msg(txt):
    global verbose
    if verbose: print >> sys.stderr, txt


        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('outfile', type=str, help="name of excel file to create")
    parser.add_argument('csv_files', type=str,nargs='+', help="list of csv file to read into the ecel file")
    parser.add_argument('--silent', '-s', action='store_true', default=False, help="silent mode (being verbose is the default)")
    parser.add_argument('--nocharts', action='store_true', default=False, help="don't create the pie charts")
    parser.add_argument('--overwrite', action='store_true', default=False, help="overwrite outfile if it exists (default is to exit)")


    args = parser.parse_args()
        
    outfile = args.outfile
    csv_files  =args.csv_files
    verbose = (args.silent == False)
    nocharts = args.nocharts
    overwrite = args.overwrite

    if os.path.exists(outfile) and not overwrite:
        msg('Error, file exists <%s>' % outfile)
        sys.exit(0)


    msg('writing Excel file %s' % outfile)
    workbook = xlsxwriter.Workbook(outfile, {'strings_to_numbers': True})
    row = -1
    
    for csv_file in csv_files:
        basename = os.path.splitext(os.path.split(csv_file)[-1])[0][:31]
        msg('\twriting worksheet %s' % basename)
        worksheet = workbook.add_worksheet(basename)

        for row, line in enumerate(open(csv_file,"r")):
            line = line.decode("UTF-8")    # Added by koala Tue 13 October 2015 (Makes it so xlsx writer knows to write unicode to output)
            fields = line.split('\t')
            newrow = row+1
            worksheet.write_row(newrow, 0, fields)

        if row == -1: continue
        if nocharts: continue
        
        if basename == 'insertSizeDistrib':
            chart = workbook.add_chart({'type': 'scatter','subtype': 'straight'})
            msg('\tdrawing scatter chart %s' % basename)
        else:
            chart = workbook.add_chart({'type': 'pie'})
            msg('\tdrawing pie chart %s' % basename)


        if basename == 'insertSizeDistrib':
            header = ['Insert Size', 'Frequency', 'Ave insert size']
            worksheet.write_row('A1', header)
            categories_col = 'A'
            values_col = 'B'
            max_rows = row + 2
            worksheet.write('C2','=SUMPRODUCT(A2:A%d,B2:B%d)/SUM(B2:B%d)' % (max_rows, max_rows, max_rows))


            chart.add_series({
                'name':       basename,
                'categories': '=%s!%s2:%s%d' % (basename, categories_col, categories_col, max_rows),
                'values':     '=%s!%s2:%s%d' % (basename, values_col, values_col, max_rows),
                })

            chart.set_x_axis({'name': 'Insert size'})
            chart.set_y_axis({'name': 'Frequency'})
            
            chart.set_size({'width': 720, 'height': 576})
            chart.set_title({'name': 'Infered insert size distribution'})
            chart.set_style(10)
            worksheet.insert_chart('E1', chart)
        elif basename.startswith('abundance.databases'):
            header = ['Mapping mode', 'database', 'Percentage' , 'Number of reads mapped']
            worksheet.write_row('A1', header)
            categories_col = 'B'
            values_col = 'C'
            max_rows = row+2
            if row > 51: max_rows = 51
            
            chart.add_series({
                'name':       basename,
                'categories': '=%s!%s3:%s%d' % (basename, categories_col, categories_col, max_rows),
                'values':     '=%s!%s3:%s%d' % (basename, values_col, values_col, max_rows),
                })
            
            chart.set_size({'width': 720, 'height': 576})
            chart.set_title({'name': basename})
            chart.set_style(10)
            
            worksheet.insert_chart('G1', chart)
        elif basename.startswith('redundant.'):
            header = ['Entry name', 'Identical entries removed from homology reduced database']
            worksheet.write_row('A1', header)

            chart.set_size({'width': 720, 'height': 576})
            chart.set_title({'name': basename})
            chart.set_style(10)
        elif basename.startswith('positive.strain'):
            header = ['Database', 'Name', 'S_Abundance (%)','R_Abundance (%)', 'Size (bp)', 'Seq_count', 'Nucleotides', 'Covered_positions', 'Coverage', 'Depth', 'Reads', 'Reads_uniq', 'Edit_distance','Description', 'Superkingdom', 'SuperkingdomTax', 'Phylum', 'PhylumTax','Class', 'ClassTax','Order', 'OrderTax','Family', 'FamilyTax','Genus', 'GenusTax','Species', 'SpeciesTax','Strain','StrainTax']
            worksheet.write_row('A1', header)
            categories_col = 'B'
            values_col = 'D'
            max_rows = row+2
            if row > 51: max_rows = 51
            chart.add_series({
                'name':       basename,
                'categories': '=%s!%s2:%s%d' % (basename, categories_col, categories_col, max_rows),
                'values':     '=%s!%s2:%s%d' % (basename, values_col, values_col, max_rows),
                })
            
            chart.set_size({'width': 720, 'height': 576})
            chart.set_title({'name': basename})
            chart.set_style(10)
            
            worksheet.insert_chart('K1', chart)
        elif basename.startswith('positive.species'):
            header = ['Database', 'Name', 'S_Abundance (%)','R_Abundance (%)', 'Size (bp)', 'Seq_count', 'Nucleotides', 'Covered_positions', 'Coverage', 'Depth', 'Reads', 'Reads_uniq', 'Edit_distance','Description', 'Superkingdom', 'SuperkingdomTax', 'Phylum', 'PhylumTax','Class', 'ClassTax','Order', 'OrderTax','Family', 'FamilyTax','Genus', 'GenusTax','Species', 'SpeciesTax']
            worksheet.write_row('A1', header)
            categories_col = 'B'
            values_col = 'D'
            max_rows = row+2
            if row > 51: max_rows = 51
            chart.add_series({
                'name':       basename,
                'categories': '=%s!%s2:%s%d' % (basename, categories_col, categories_col, max_rows),
                'values':     '=%s!%s2:%s%d' % (basename, values_col, values_col, max_rows),
                })
            
            chart.set_size({'width': 720, 'height': 576})
            chart.set_title({'name': basename})
            chart.set_style(10)
            
            worksheet.insert_chart('K1', chart)
        elif basename.startswith('negative.strain'):
            header = ['Database', 'Name', 'S_Abundance (%)','R_Abundance (%)', 'Size (bp)', 'Seq_count', 'Nucleotides', 'Covered_positions', 'Coverage', 'Depth', 'Reads', 'Reads_uniq', 'Edit_distance','Description', 'Superkingdom', 'SuperkingdomTax', 'Phylum', 'PhylumTax','Class', 'ClassTax','Order', 'OrderTax','Family', 'FamilyTax','Genus', 'GenusTax','Species', 'SpeciesTax','Strain','StrainTax']
            worksheet.write_row('A1', header)
            categories_col = 'B'
            values_col = 'D'
            max_rows = row+2
            if row > 51: max_rows = 51
            chart.add_series({
                'name':       basename,
                'categories': '=%s!%s2:%s%d' % (basename, categories_col, categories_col, max_rows),
                'values':     '=%s!%s2:%s%d' % (basename, values_col, values_col, max_rows),
                })
            
            chart.set_size({'width': 720, 'height': 576})
            chart.set_title({'name': basename})
            chart.set_style(10)
            
            worksheet.insert_chart('K1', chart)
        elif basename.startswith('negative.species'):
            header = ['Database', 'Name', 'S_Abundance', 'R_Abundance (%)', 'Size (bp)', 'Seq_count', 'Nucleotides', 'Covered_positions', 'Coverage', 'Depth', 'Reads', 'Reads_uniq', 'Edit_distance','Description', 'Superkingdom', 'SuperkingdomTax', 'Phylum', 'PhylumTax','Class', 'ClassTax','Order', 'OrderTax','Family', 'FamilyTax','Genus', 'GenusTax','Species', 'SpeciesTax']
            worksheet.write_row('A1', header)
            categories_col = 'B'
            values_col = 'D'
            max_rows = row+2
            if row > 51: max_rows = 51
            chart.add_series({
                'name':       basename,
                'categories': '=%s!%s2:%s%d' % (basename, categories_col, categories_col, max_rows),
                'values':     '=%s!%s2:%s%d' % (basename, values_col, values_col, max_rows),
                })
            
            chart.set_size({'width': 720, 'height': 576})
            chart.set_title({'name': basename})
            chart.set_style(10)
            
            worksheet.insert_chart('K1', chart)

    workbook.close()
