#! /usr/bin/python -u
__author__  = "Marco Mariotti"
__email__   = "mmariotti@bwh.harvard.edu"
__licence__ = "GPLv3"
__version__ = "1.0"
from string import *
import sys
import traceback
sys.path.append('/home/mmariotti/scripts')
from MMlib import *
from tree_classes import syntheny_view
from ete2 import Tree, TreeStyle, NodeStyle, faces
import random

help_msg="""Program to show a graphical representation based on ETE2 of the syntheny around some genes of interest. Each gene is represented as a colored arrow.
Usage:  $ syntheny_view.py  -i genes.gff -a annotation.gff -f homology.tsv  [options]  

##### Input options (all required):
-i      gff annotation file of genes of interest 
-a      global gff annotation file. Genes overlapping with those of interest will be removed
-f      homology tab separated file (i.e. lines like "geneId -tab- familyId"); to obtain one, run a all-against-all blastp and then run blast_homology_clusters.py 

## processing input options
-if     function to apply to each line of the -i file to determine the gene id. Default: first word of last field -- i.e. -if "x.split('\t').split()[0]" 
-af     same as -if, but for the -a file. 
-at     tag (element type, third tab field) of the lines that will be kept from the -a file. Default: "gene";  use "*" to keep every element

## gene windows options
-l      length of nts on either side of each gene of interest to be displayed
-n      number of genes on either side of each gene to be displayed (overrides -l)

## colors options
-c      file with available colors. File with a single word per line, each identified as color in ete2/Qt4 (e.g. "#4567A1", "turquoise"). These colors are used as many times as necessary to draw all families
-cr     color rotation scheme; defines how to use the available colors. Arg must be a number between 0 and 3, referring cumulatively to [0:fill, 1:outline, 2:box_bkg, 3:box_outline]. Example, with a value of 1 (default), colors are used for the fill and the outline only, while the boxes are kept always transparent. Use a higher -cr when in need to represent many families
-rc     randomize the colors in the -c file
-cs     color for singlets, i.e. genes without a family assigned, or those being the only representatives for their family in the current run. Arg accepted: like -ci. Grey-like colors are suggested. Default: "gainsboro,darkgrey"
-cf     color-per-family file. First word of each line must be a familyId, present in the input homology tsv file (or, a geneId). The next words (at least one must be present) are taken respectively as fill color, outline color, box background color, box outline color. This option overrides the colors that would be chosen with options -c and -cr. You can provide even just a subset of all families, with the rest being drawn according to -c and -cr
-ci     color of genes of interest. You can provide a single value (fill) or comma separated values (max 4), which are interpreted as fill, outline, box background and box outline. This option has priority on all other color options. You should use colors that are not among the available colors provided with -c. Default: "white,royalblue". 

## other graphical attributes
-w       width in pixels of each gene arrow (default:120)
-fs      font size used  (default:8)
-m       minimal output; information in each gene arrow is reduced. Also, -w is set to 30

## modifications of the gene cluster representation
-rs      remove all singlets before displaying (see -cs for singlet definition)
-rf      remove families; provide a file one familyId (geneId also accepted) per line. These will be ignored from the annotation file
-ks      keep strand; by default, if the gene of interest is on the - strand, everything is flipped. Use this option to keep it as it is
-dm      don't merge; by default, if the resulting clusters share some genes, they are joined before displaying. Use this option to avoid this, so to have always a gene cluster per input gene

## output text files; provide an argument that will be the output file
-of      family output; to get the list and composition of families in the current view. One line per family, like "family geneId1 geneId2 etc" (with tabs)
-oc      clusters output; to get the gene clusters in the current view. One line per cluster, like "geneId1 geneId2 etc" (with tabs)
-ocf     color per family output; one line per family in current view, like "family color_bkg color_outline box_bkg box_outline" (with tabs), using "None" for undefined colors 
-ocg     color per gene output; one line per gene in current view, like "geneId color_bkg color_outline box_bkg box_outline" (with tabs), using "None" for undefined colors  

## miscellaneous
-out            instead of opening the interactive ete2 environment (default), create this output file (pdf or png extensions are accepted in the arg)
-temp           temporary folder; a subfolder is created here and deleted upon exiting
-legend         suppress normal input and output. Provide a file with lines like "left_id -tab- color_bkg color_outline box_bkg box_outline -tab- arrow_text -tab- desc_text". This will build a legend-like representation with one arrow per line colored as specified (you can use None as color), some text inside the arrow, and some next to it. Strings "\n" will be interpreted as newline characters.
-print_opt      print currently active options
-v              verbose; prints a lot of stuff
-h OR --help    print this help and exit"""

command_line_synonyms={}

def_opt= { 'temp':'/home/mmariotti/temp', 
'i':'',   'a':'',  'f':'',
'if':None,    'af':"x.split('\t')[-1].split(';Name=')[-1].split(';')[0]",     'at':'gene',
'c':'/home/mmariotti/libraries/available_colors.tab', 'cr':1,  'cf':0,  'ci': 'white,royalblue', 'cs':'gainsboro,darkgrey', 'rc':0,
'l':10000, 'n':0,
'fs':8, 'w':120,
'm':0, 
'rs':0, 'ks':0, 'dm':0, 'rf':0,
'of':0, 'oc':0, 'ocf':0, 'ocg':0,
'v':0, 'out':0,
'legend':0,
}

#########################################################
###### start main program function

class gene_cluster(list):
  """ Simple class to add  a single attribute to a list of genes, which is: the central gene of interest used to populate this list (.ref_gene)"""
  def link_to_gene(self, g):     self.ref_gene=g
  
class notracebackException(Exception):
  """ when raising one of these, the python error message will be much less verbose than the standard traceback """

def main(args={}):
#########################################################
############ loading options
  global opt
  if not args: opt=command_line(def_opt, help_msg, 'iaf', synonyms=command_line_synonyms, strict=1)
  else:  opt=args
  set_MMlib_var('opt', opt)
  global temp_folder; temp_folder=Folder(random_folder(opt['temp'])); test_writeable_folder(temp_folder, 'temp_folder'); set_MMlib_var('temp_folder', temp_folder)

  write("#=============------ "); write(" show_syntheny.py ", how='reverse,bright'); write(" -------=============#", 1)
  ## basic graphics options
  face_width= opt['w'];   font_size= opt['fs']
  ## defining a function that, given a gene, decide what is printed in the face"""
  if not opt['m']: 
    def get_text(g): 
      if g.id in geneid2family:      return geneid2family[g.id]+':'+g.id
      else:                          return '-'+':'+g.id
  else:
    face_width=30
    def get_text(g):                 
      if g.id in geneid2family:      return geneid2family[g.id]
      else:                          return ''
  tree_style=TreeStyle(); tree_style.show_leaf_name=False; tree_style.scale=1; tree_style.show_scale=False
  node_style=NodeStyle(); node_style["size"] = 0 #; node_style["fgcolor"] = "darkred"
  node_style_grey=NodeStyle(); node_style_grey["size"] = 0; node_style_grey["bgcolor"] = "lightgrey"
  tree=Tree(); tree.dist=0.0;  tree.set_style(node_style)

  ############################## legend mode only: start
  if opt['legend']: ### totally different program in this case
    for line_index, line in enumerate(open(opt['legend'])):
      try:
        bkg_color='white' if line_index%2 else 'lightgrey' 
        splt=line.strip().split('\t')
        if not splt: continue
        leaf=tree.add_child(name='', dist=0.0);   
        leaf.set_style(node_style) if line_index%2 else leaf.set_style(node_style_grey)
        g=gene(strand='+'); g.color, g.color_outline, g.color_box_bkg, g.color_box_line=[x if x!='None' else None for x in splt[1].split()]; 
        g.text = replace(splt[2], '\\n', '\n')
        title_face=faces.TextFace(splt[0], fsize=font_size); title_face.margin_left=5;  leaf.add_face(title_face, 0, 'branch-right' ) #title left
        arrow_face=syntheny_view([g],   printed={'boundaries': 0, 'text':1, 'id':0},  pen_size=4, font_size=font_size, width=face_width)[0]; leaf.add_face(arrow_face, 1, 'branch-right' )
        for desc_line_index, desc_line in enumerate(splt[3].split('\\n')):
          desc_face=faces.TextFace(desc_line, fsize=font_size);  desc_face.background.color = bkg_color; leaf.add_face(desc_face,  2, 'branch-right' ) #desc_face.margin_down=3; desc_face.margin_up=3;
      except: printerr('-legend ERROR parsing this line: |n{0}'.format(line), 1); raise
    write("Legend mode: {0} entries found. ".format(len(tree)), 1)
    if opt['out']:    write('--> writing output file: {0}'.format(opt['out']), 1);       tree.render(opt['out'], tree_style=tree_style)
    else:             write('-- opening interactive ETE2 environment (PyQt4) -- ', 1);  tree.show(tree_style=tree_style)
    sys.exit()
  ############################## legend mode only: over

  #### checking input
  input_gff_file=opt['i'];      check_file_presence(input_gff_file, 'input_gff_file', notracebackException )
  annotation_gff_file=opt['a']; check_file_presence(annotation_gff_file, 'annotation_gff_file', notracebackException )
  homology_file=opt['f'];       check_file_presence(homology_file, 'homology_file', notracebackException )
  
  # printing for pretty out
  write('# Input gff file=      {0:<30} (genes of interest)'.format(input_gff_file), 1)
  write('# Annotation gff file= {0:<30} (all genes)'.format(annotation_gff_file), 1)
  write('# Homology tsv file=   {0:<30} (gene families)'.format(homology_file), 1)
  non_def_options_str=join([ '# -{0}  {1}\n'.format(k, opt[k]) for k in opt  if k in def_opt and def_opt[k] != opt[k] and not k in 'iaf' ], '')
  if non_def_options_str:  write('### Non-default options:\n'+non_def_options_str)
  write('', 1)

  # checking output options
  for x in ['of', 'oc', 'ocf', 'ocg']:
    if opt[x] and opt[x]==1: raise notracebackException, "ERROR option -{0} must be provided with an argument (which will be used as output file)!"

  #######
  ### processing options controlling colors
  colors_already_taken={}  # useful for later, when we compute available_colors
  color_genes_of_interest=[None, None, None, None]
  if opt['ci']: 
    for index, color in enumerate( opt['ci'].split(',') ):  
      if color=='None': color=None
      color_genes_of_interest[index]=color
    colors_already_taken[ join( map(str, color_genes_of_interest), ',') ]=1

  color_singlets=[None, None, None, None]
  if opt['cs']: 
    for index, color in enumerate( opt['cs'].split(',') ):  
      if color=='None': color=None
      color_singlets[index]=color
    colors_already_taken[ join( map(str, color_singlets), ',') ]=1

  fam2color={}    ## each color is a list [fill, outline, box_bkg, box_outline]  if not defined, it's None
  if opt['cf']:
    ## load color-family file
    for line in open( opt['cf'] ):
      splt=line.strip().split()
      if splt:  #skipping empty lines
        fam=splt[0];  the_colors=[None, None, None, None]
        for index, item in enumerate(splt[1:]): 
          if item=='None': item=None
          the_colors[index]=item
        fam2color[fam] = the_colors
        colors_already_taken[ join( map(str, the_colors), ',') ]=1

  color_file=opt['c'];          check_file_presence(color_file, 'color_file', notracebackException )
  color_scheme=opt['cr'];       
  if not color_scheme in [0, 1, 2, 3]: raise notracebackException, "ERROR invalid color scheme provided with option -cr ! see -help"
  individual_colors=[line.strip() for line in open(color_file) if line.strip()]; 
  if     color_scheme==0:  available_colors = [[a,None,None,None] for a in individual_colors if not a+',None,None,None' in colors_already_taken]
  elif   color_scheme==1:  available_colors = [[b,   a,None,None] for a in individual_colors for b in individual_colors if not b+','+a+',None,None' in colors_already_taken]
  elif   color_scheme==2:  available_colors = [[c,   b,   a,None] for a in individual_colors for b in individual_colors for c in individual_colors if not (a==b==c) and not c+','+b+','+a+',None' in colors_already_taken]
  elif   color_scheme==3:  available_colors = [[d,   c,   b,   a] for a in individual_colors for b in individual_colors for c in individual_colors for d in individual_colors if not (b==c==d) and not d+','+c+','+b+','+a in colors_already_taken]
  #write('available colors: {0}'.format(len(available_colors)), 1)

  if opt['rc']:   random.shuffle(available_colors)

  ######
  ## loading gff input files    # genes of interest
  input_get_id_function=None; 
  if opt['if']: input_get_id_function=eval('lambda x:'+opt['if'])
  write('Loading genes of interest from {0:<30} ... '.format(input_gff_file)) 
  genes_of_interest=load_all_genes(input_gff_file, tag='*', get_id=input_get_id_function, is_sorted=True)
  for g_index, g in enumerate(genes_of_interest): g.is_of_interest=g_index+1   ### keeping this as a number so later we can sort output in the same order as input
  write('done. Genes: {0}'.format(len(genes_of_interest)), 1)
                                # gene in global annotation
  annotation_get_id_function=None; 
  if opt['af']: annotation_get_id_function=eval('lambda x:'+opt['af'])
  annotation_tag=opt['at']
  write('Loading annotated genes from   {0:<30} ... '.format(annotation_gff_file)) 
  annotated_genes=load_all_genes(annotation_gff_file, tag=annotation_tag, get_id=annotation_get_id_function)
  for a in annotated_genes: a.is_of_interest=False
  write('done. Genes: {0}'.format(len(annotated_genes)), 1)
  ######

  ## load homology file
  geneid2family={}; families_dict={}
  write('Loading homology families from {0:<30} ... '.format(homology_file)) 
  for line in open(homology_file):
    splt=line.strip().split('\t')
    if splt:  geneid, family = splt; geneid2family[geneid]=family; families_dict[family]=0
  write('done.', 1)

  ## print some stats
  for g in annotated_genes: 
    if g.id in geneid2family:  families_dict[geneid2family[g.id]]+=1
  n_fam_represented=0; n_genes_with_family=0
  for fam in families_dict: 
    if families_dict[fam]>0: n_fam_represented+=1; n_genes_with_family+=families_dict[fam]
  write('N of families: {0} ; {1} families have 1 or more gene(s) found in annotation.\nA total of {2} genes have a family assigned.\n'.format(len(families_dict),n_fam_represented, n_genes_with_family ), 1)
  del families_dict;  #saving memory (almost a joke)
  family2genes_displayed={}      ### later we'll modify geneid2family to avoid displaying useless families

  ## families or genes in the annotation to be ignored
  fams_to_ignore={}
  if opt['rf']: 
    check_file_presence(opt['rf'], '-rf file')
    for line in open(opt['rf']):
      strp=line.strip()
      if strp:  fams_to_ignore[strp]=True

  ##############################  start doing things!
  ## finding overlaps
  def scoring_function_for_overlaps(g):    return int (g.is_of_interest)* 10000000 + g.length()
  removed_overlapping_genes=[]
  non_red_genes = remove_overlapping_gene_clusters( genes_of_interest + annotated_genes,  scoring=scoring_function_for_overlaps, phase=True, strand=True, out_removed_genes=removed_overlapping_genes, remember_overlaps=True )
  ### getting all discarded -> kept  relationship, and back
  for g in removed_overlapping_genes: 
    if not hasattr( g.overlapping, 'discarded'): g.overlapping.discarded=[]
    g.overlapping.discarded.append( g )
  for g in genes_of_interest: 
    if hasattr( g, 'discarded'): #len(g.discarded)>1: 
      for d in g.discarded:   write(' Gene: {0:^25} removed overlapping gene: {1}'.format(g.id, d.id), 1)

  non_red_genes.sort(  cmp=order_genes_for_chr_pos  )   #sorting again... not optimized but easy
  ######

  ##############################
  ## building gene clusters to be displayed
  gene_clusters=[]   # list of lists of genes; populating this while parsing the sorted list of genes and looking for the genes of interest.
  max_distance = opt['l']
  index=0
  while index < len(non_red_genes):
    if non_red_genes[index].is_of_interest:
      g=non_red_genes[index]
      verbose('*** Cluster of {0}'.format(g.id), 1)
      gc=gene_cluster();  gc.append(g); gc.link_to_gene(g)

      ## parsing back CAREFUL ASSUMING THERE ARE NO NESTED STRUCTURES WITH EXONS
      other_index=index-1
      while other_index >= 0 and non_red_genes[other_index].chromosome == g.chromosome and  \
        (  (not opt['n'] and abs( g.boundaries()[0] - non_red_genes[other_index].boundaries()[1] ) <= max_distance )  or \
           (    opt['n'] and len(gc)-1 <= opt['n']    )    ):
        if not  non_red_genes[other_index].id in fams_to_ignore and not \
           (non_red_genes[other_index].id in geneid2family and geneid2family[non_red_genes[other_index].id] in fams_to_ignore):
          gc.insert(0, non_red_genes[other_index])
        other_index-=1
      n_genes_added_back=len(gc)-1

      #parsing forward                                                                                               
      other_index=index+1
      while other_index < len(non_red_genes) and non_red_genes[other_index].chromosome == g.chromosome and \
        (   (not opt['n'] and abs( non_red_genes[other_index].boundaries()[0] - g.boundaries()[1] ) <= max_distance )  or \
            (    opt['n'] and len(gc)-1-n_genes_added_back <= opt['n']    )    ):
        if not  non_red_genes[other_index].id in fams_to_ignore and not \
           (non_red_genes[other_index].id in geneid2family and geneid2family[non_red_genes[other_index].id] in fams_to_ignore):            
          gc.append(non_red_genes[other_index])
        other_index+=1

      for i in gc:       verbose( i.gff(), 1)
      gene_clusters.append(gc) 

    index+=1

  ## populating family2genes_displayed to compress family output
  for gc in gene_clusters:
    for g in gc: 
      if g.id in geneid2family: 
        fam=geneid2family[g.id]
        if not fam in family2genes_displayed: family2genes_displayed[fam]={}
        family2genes_displayed[fam][g.id]=True

  if opt['rs']:
    ## removing singlets
    n_singlets_removed=0
    for gc in gene_clusters:
      len_gc=len(gc)
      for i in range(len_gc): 
        g_index= len_gc-i-1  #parsing in reverse order to make .pop() work
        g= gc[g_index]
        if  not g.is_of_interest and ( not g.id in geneid2family or len(family2genes_displayed[ geneid2family[g.id] ])==1 ):  
          gc.pop(g_index); n_singlets_removed+=1;  del family2genes_displayed[ geneid2family[g.id] ]
    if n_singlets_removed: write('Option -rs: {0} singlets were removed! '.format(n_singlets_removed), 1)

  #### merging clusters that share at least one gene
  if not opt['dm']:
    ## since they are sorted, a cluster can share genes only with its previous or next cluster. Also, if we scan forward, we just check if the last gene in a cluster is contained in the next one
    gc_index=0
    while gc_index+1<len(gene_clusters):   #+1 since, if it's the last one, it's not interesting
      current_gc=gene_clusters[gc_index];    next_gc=gene_clusters[gc_index+1]
      merged=False
      if current_gc[0].chromosome == next_gc[0].chromosome: # and current_gc[-1] in next_gc:     --> in practice this is what we check. but let's do it more efficiently
        try:     
          index_pos= next_gc.index(current_gc[-1])  # this cause an exception if not there
          merged=True             #### Yes we're officially merging
          write('Merging the surrounds of {0:>25} and {1:<25}'.format(current_gc.ref_gene.id,  next_gc.ref_gene.id), 1)
          for gc in next_gc[index_pos+1:]:   current_gc.append( gc )
          possible_ref_genes=[]
          for g_index, g in enumerate(current_gc):
            if g.is_of_interest: g.g_index=g_index; possible_ref_genes.append(g)
          middle_point = (len(current_gc)-1)/2.0
          best_ref_gene = min (possible_ref_genes, key= lambda x:abs(x.g_index-middle_point))
          current_gc.link_to_gene(best_ref_gene)
          gene_clusters.pop(gc_index+1)   #removing next_gc
        except ValueError: pass
      if not merged:    gc_index+=1
       
  geneid2color={}
  ### parsing each single gene in each cluster. deciding COLORS
  for gc in gene_clusters:  
    for g in gc: 
      if g.id in fam2color:          geneid2color[g.id]=fam2color[g.id]    # color was specified in -cf using geneId      
      elif g.id in geneid2family:   
        ## this belongs to a family
        fam = geneid2family[g.id]  
        if len(family2genes_displayed[fam]) == 1:   
          geneid2color[g.id]= color_singlets   ## singlet being only representative for its family
        else:   
          if not fam in fam2color:     # not yet assigned to this family
            try:                fam2color[fam]=available_colors.pop(0)
            except IndexError:  raise notracebackException, "ERROR not enough colors are available to display this! Increase the number of colors in the -c file or change the color scheme with -cr"
          geneid2color[g.id]= fam2color[fam]
      else:    geneid2color[g.id]=color_singlets     ## singlet that does not belong to any family
      if g.is_of_interest:
        geneid2color[g.id]= list(geneid2color[g.id])  ## copying list or otherwise we modify in place the color
        for index, color in enumerate(color_genes_of_interest):
          if not color is None:  geneid2color[g.id][index]=color
    #write( g.id+' '+str(geneid2color[g.id]) +'    '+ str(geneid2color), 1, how='green')
    #write('---', 1, how='reverse')

  ### now sorting gene_clusters so they are in the same order as the input file
  gene_clusters.sort(key=lambda x:x.ref_gene.is_of_interest)

  if opt['of']:
    ### producing an output file with a line for each family
    fh=open(opt['of'], 'w')
    for fam in family2genes_displayed: print >> fh, fam+'\t'+join( [gid for gid in family2genes_displayed[fam]], '\t')
    fh.close()
  if opt['oc']:
    ### producing an output file with a line for each gene cluster
    fh=open(opt['oc'], 'w')
    for gc in gene_clusters: print >> fh, join( [g.id for g in gc], '\t')
    fh.close()
  if opt['ocf']:
    fh=open(opt['ocf'], 'w')
    for fam in fam2color:
      if not fam in family2genes_displayed: continue #this is to skip cases in which -cf was provided with geneId instead of fam
      print >>fh,  fam+'\t'+join(map(str, fam2color[fam]), '\t')
    fh.close()
  if opt['ocg']:
    fh=open(opt['ocg'], 'w')
    for geneid in geneid2color: print >>fh,  geneid+'\t'+join(map(str, geneid2color[geneid]), '\t')
    fh.close()


  ### preparing ete2 objects
  max_n_genes_in_a_cluster= max ([len(gc) for gc in gene_clusters])
  for gc in gene_clusters:
    name_displayed= gc[0].chromosome + ' : ' +gc.ref_gene.id
    name_displayed+=join([ '\n'+ ' & '+g.id    for g in gc if g.is_of_interest and g != gc.ref_gene], '') #adding a line for other genes of interested merged in this cluster
    leaf=tree.add_child(name=name_displayed, dist=10);   leaf.set_style(node_style)
    leaf_name_face= faces.TextFace(leaf.name, fsize=font_size);     leaf_name_face.margin_left=5;  leaf_name_face.margin_right=1; 
    leaf.add_face( leaf_name_face, 0, 'aligned' )
    for g in gc:  
      g.color, g.color_outline, g.color_box_bkg, g.color_box_line  = geneid2color[g.id]
      g.text = get_text(g)

    reverse_syntheny_view =  not opt['ks'] and gc.ref_gene.strand=='-'
    # modifying gc inplace to add whitespacers to keep it sortof centered
    while len(gc)<max_n_genes_in_a_cluster:     
      if (len(gc) + int(reverse_syntheny_view))%2 :   gc.append('')
      else:           gc.insert(0, '')
    face_list=     syntheny_view(  gc,   printed={'boundaries': int(not opt['m']), 'text':1, 'id':0},  pen_size=4, font_size=font_size, width=face_width, reverse=reverse_syntheny_view)
    for col_index, the_face in enumerate(face_list):   leaf.add_face( the_face, col_index+1, 'aligned' )
  
  if opt['out']: 
    write('--> writing output file: {0}'.format(opt['out']), 1)
    tree.render(opt['out'], tree_style=tree_style)
  else:          
    write('-- opening interactive ETE2 environment (PyQt4) --- ', 1)
    tree.show(tree_style=tree_style)

  write('#====----      execution finished, exiting         -----====#', 1)

  ###############



#######################################################################################################################################

def close_program():
  if 'temp_folder' in globals() and is_directory(temp_folder):
    bbash('rm -r '+temp_folder)
  try:
    if get_MMlib_var('printed_rchar'): 
      printerr('\r'+printed_rchar*' ' ) #flushing service msg space       
  except:
    pass

  if sys.exc_info()[0]: #an exception was raised. let's print the information using printerr, which puts it in the logfile as well, if any.   
    if issubclass(sys.exc_info()[0], notracebackException):      printerr( sys.exc_info()[1], 1)
    elif issubclass(sys.exc_info()[0], SystemExit):      pass
    else:                                                        printerr('ERROR '+ traceback.format_exc( sys.exc_info()[2]) , 1)

  if 'log_file' in globals(): log_file.close()


if __name__ == "__main__":
  try:
    main()
    close_program()  
  except Exception:
    close_program()
