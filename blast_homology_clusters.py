#! /usr/bin/python -u
from string import *
import sys
from commands import *
#sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/home/mmariotti/scripts')
#sys.path.append('/home/mmariotti/software/selenoprofiles')
from MMlib import *
from annotate_with_tblastn import parse_blast_tab

help_msg=""" Parse tabular output (m8) of blast and builds homology families. If two proteins have a blast hit satisfying the filter, they are joined in a family. Two proteins can be in the same family without having a blast hit linking them if there's a protein with blast hits to both (single link clustering). 
Normal output is human readable (although may be very long). For usage with other programs (e.g. syntheny_view.py) see option -A

Usage: $  blast_homology_clusters.py  blast_output.tab  [options]  > output

-e      blast evalue
-b      blast output is default blastall, not tabular

-n      require at least N members in a cluster to output a family
-m      in normal output (no -A) defined max examples shown for cluster

-A      tab output in families, as for the adhore program. Example of a line:   protein1-tab-F1
-add    provide fasta file to add proteins with no hits as single member families (if -A is active) or for stats in normal output

-s      do not count links involving same species
-sf     species function, in lambda style. Every protein name is evaluated in this way to extract species name 

### Options:
-print_opt      print currently active options
-h OR --help    print this help and exit"""

command_line_synonyms={}

def_opt= {#'temp':'/users/rg/mmariotti/temp', 
'i':0, 
'e':'1e-10', 'b':0, 
'A':0, 'add':0, 
'm':5, 
  's':1, 'sf': "x.split('.')[0]",
'v':0, 'n':0,
}


#########################################################
###### start main program function

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'io', synonyms=command_line_synonyms )
  else:  opt=args
  set_MMlib_var('opt', opt)
  #global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)
  #global split_folder;    split_folder=Folder(opt['temp']);               test_writeable_folder(split_folder); set_MMlib_var('split_folder', split_folder) 
  #checking input
  input_file=opt['i'];   check_file_presence(input_file, 'input_file')
  all_titles_short2long= {} 
  if opt['add']: all_titles_short2long=   dict( (t.split()[0], t)  for t,s in parse_fasta(opt['add']) )
  max_examples=opt['m']
  evalue_threshold=e_v(opt['e'])
  species_function= eval('lambda x:'+opt['sf'])

  cluster_index2geneids={}; cluster_index=1
  geneid2cluster_index={}

  for bhit in parse_blast_tab(input_file):
    if not  bhit.evalue < evalue_threshold: continue

    id_left, id_right=     bhit.chromosome, bhit.query.chromosome
    if opt['s']:           species_left, species_right=  species_function( id_left ) , species_function( id_right )
#    species_left, species_right=  id_left.split('.')[0], id_right.split('.')[0]
    if   id_left != id_right    and (not opt['s'] or species_left != species_right):
      if   (  not id_left in geneid2cluster_index )  and  ( not id_right in geneid2cluster_index ):  #new cluster
        cluster_index2geneids [cluster_index] = [id_left, id_right]
        geneid2cluster_index[id_left]=cluster_index;              geneid2cluster_index[id_right]=cluster_index
        cluster_index+=1
#        print '1 creating with ids: ',id_left, id_right, ' the cluster ', cluster_index-1
      elif (  id_left in geneid2cluster_index     )  and  ( not id_right in geneid2cluster_index ):
        cluster_index2geneids [ geneid2cluster_index[id_left] ].append( id_right )
        geneid2cluster_index[id_right]=   geneid2cluster_index[id_left]
#        print '2 moving id: '+id_right+' to cluster ', cluster_index
      elif ( not id_left in geneid2cluster_index  )  and  ( id_right in geneid2cluster_index ):
        cluster_index2geneids [ geneid2cluster_index[id_right] ].append( id_left )
        geneid2cluster_index[id_left]=    geneid2cluster_index[id_right]
#        print '3 moving id: '+id_left+' to cluster ', cluster_index
      elif geneid2cluster_index[id_left] != geneid2cluster_index[id_right]    : #not id_left in geneid2cluster_index  and  not id_right in geneid2cluster_index   is implicit
        #putting those in cluster of id_right into the cluster in id_left, unless the reverse is more efficient
        if len( cluster_index2geneids [ geneid2cluster_index[id_right] ] )  > len( cluster_index2geneids [ geneid2cluster_index[id_left] ] ): id_left, id_right = id_right, id_left
#        print '4 LEFT:', id_left,  'c:', geneid2cluster_index[id_left],  'RIGHT:',  id_right, 'c:', geneid2cluster_index[id_right]

        cluster_index_to_remove = geneid2cluster_index[id_right]
        for gid in cluster_index2geneids [ geneid2cluster_index[id_right] ]:
#          print '4 moving id: '+gid, 'from cluster', cluster_index_to_remove, ' to cluster ', geneid2cluster_index[id_left]
          geneid2cluster_index[ gid ] =  geneid2cluster_index[id_left]
        cluster_index2geneids[  geneid2cluster_index[id_left]  ].extend(   cluster_index2geneids[  cluster_index_to_remove  ]   )
#        print '4 removing cluster', cluster_index_to_remove
        del cluster_index2geneids[  cluster_index_to_remove  ]


  ordered_cluster_indexes=sorted( cluster_index2geneids.keys(), key=lambda x:len(cluster_index2geneids[x]), reverse=True    )  #largest clusters first

  if opt['n']:   
    for fam_index,  cluster_index in enumerate(ordered_cluster_indexes[::-1]):
      if len(cluster_index2geneids[cluster_index]) < opt['n']: 
        del cluster_index2geneids[cluster_index]
        ordered_cluster_indexes.pop( -1 )

  ###### now it's time to output
  if opt['A']:
    for fam_index, cluster_index in enumerate(ordered_cluster_indexes):
      #printerr('', 1)
      for gid in cluster_index2geneids[cluster_index]:
        write('{0}\tF{1}'.format(gid, fam_index+1), 1)
  else:
    write('N clusters: {0}'.format(len(cluster_index2geneids.keys())), 1)
    for fam_index, cluster_index in enumerate(ordered_cluster_indexes):
      write('Cluster F{1} has {0:^5} elements:'.format(len( cluster_index2geneids[cluster_index]), fam_index+1), 1, how='red')
      for short_t in cluster_index2geneids[cluster_index]  [:max_examples]:
        long_t=''
        if opt['add']: 
          try:    long_t= all_titles_short2long[short_t].split('#')[1] [:100]
          except: pass
        write(' '+short_t+' '+long_t , 1)


#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass

  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
    raise 
