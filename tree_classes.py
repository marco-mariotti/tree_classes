#! /usr/bin/python -u
import sys
sys.path.insert(0, "/users/rg/mmariotti/libraries/")
sys.path.append('/users/rg/mmariotti/scripts')
from MMlib import nogap
from string import join
try:  
  from PyQt4 import QtCore, QtGui
  from ete2 import faces
except: 
  printerr('ERROR PyQt4 and ete2 must be installed to use this program!', 1)
  raise 

global warnings_reduced_text; warnings_reduced_text={} #keys: categories
def shrink_font_to_box (simpleTextItem, width=-1, height=-1, category='unknown'):
  """Function to shrink a piece of text in a rectangle of a given width and height; its initial font size value is reduced until it fits. The global variable warning_reduced_text keeps track of every time this function is run succesfully. Note: after calling this function you set the positions of the simpleTextItem object according to its renewed boundingRect """
  psize=99 #not used, this is just for the first time last condition of while is checked
  reducing_of_points=0
  while (width>0 and simpleTextItem.boundingRect().width() > width) or (height>0 and simpleTextItem.boundingRect().height()> height ) and psize>1:
    font=simpleTextItem.font()
    psize= font.pointSize()
    font.setPointSize(psize-1)
    simpleTextItem.setFont(font)
    reducing_of_points+=1
  if reducing_of_points:
    if not warnings_reduced_text.has_key(category): warnings_reduced_text[category]=0
    warnings_reduced_text[category]+=1

class GapFace(faces.Face):
  """ Simple face that can be used as place holder, e.g. to indicate gaps in the syntheny view. A single horizontal dotted line is drawn at the center. Width and height can be used so that this can be aligned to other faces.""" 
  def __init__(self, width=150, height=40, color='gray', pen_size=1):
    faces.Face.__init__(self)
    self.type = "item"
    self.item = None
    self.width=width
    self.height=height
    self.color=color
    self.pen_size=pen_size
  def _width(self):        return self.item.rect().width()
  def _height(self):       return self.item.rect().height()
  def update_items(self):
      self.item = QtGui.QGraphicsRectItem(0, 0, self.width, self.height) # backbone object, no color 
      self.item.setPen(QtGui.QPen(QtCore.Qt.NoPen)) #no black border      
      if not self.color is None:
        outline_pen = QtGui.QPen(    QtGui.QColor(self.color),  self.pen_size, QtCore.Qt.DotLine)
        gap_line =  QtGui.QGraphicsLineItem( 0, self.height/2.0,  self.width,  self.height/2.0,            self.item)
        gap_line.setPen(outline_pen)

class ContigTerminusFace(faces.Face):
  """ Simple face that can be used as place holder, to indicate the end of this contig (or the start)  in the syntheny view. Width and height can be used so that this can be aligned to other faces. Init as: cef= ContigEndFace(width=100, height=40, thickness=5, arm_width=20, arm_height=30, arcolor='gray', mode='start')    OR: mode='end'
  <              width           >
  ################################ ^            mode=start                   mode=end
  #       <        > arm_width   #                                    #   
  #     ^ ---------| ^           # h        x1,y4             x3,y4   #   x3,y4               x1,y4
  #       | -------| v thickness # e              x2,y3       x3,y3   #   x3,y3         x2,y3
  # arm_  | |                    # i                                  #
  # height| |                    # g                                  #
  #       | -------|             # h              x2,y2       x3,y2   #   x3,y2         x2,y2
  #     v ---------|             # t        x1,y1             x3,y1   #   x3,y1               x1,y1
  #       < > thickness          #                                    #
  ################################ v
""" 
  def __init__(self, width=150, height=40, thickness=5, arm_width=20, arm_height=30,color='darkgray', mode='start'):
    faces.Face.__init__(self)
    self.type = "item"
    self.item = None
    if arm_width  > width :  arm_width=  width
    if arm_height > height : arm_height= height
    thickness= min([thickness, width, height])
    self.width=width;     self.height=height; 
    self.color=color; self.mode= mode
    self.thickness=thickness;   self.arm_width=arm_width;  self.arm_height=arm_height

  def _width(self):        return self.item.rect().width()
  def _height(self):       return self.item.rect().height()
  def update_items(self):
      self.item = QtGui.QGraphicsRectItem(0, 0, self.width, self.height) # backbone object, no color 
      self.item.setPen(QtGui.QPen(QtCore.Qt.NoPen)) #no black border      
      fill_brush=  QtGui.QBrush(  QtGui.QColor(self.color)  )
      #pen_size=1
      #outline_pen = QtGui.QPen(    QtGui.QColor(self.color),  pen_size, QtCore.Qt.DotLine)
      if self.mode=='start':        x1 = (self.width - self.arm_width) / 2.0;       x2 = x1 + self.thickness;       x3 = x1 + self.arm_width
      elif self.mode=='end':        x3 = (self.width - self.arm_width) / 2.0;       x1 = x3 + self.arm_width;       x2 = x1 - self.thickness
      y4 = (self.height - self.arm_height) / 2.0;     y3 = y4 + self.thickness;       y1 = y4 + self.arm_height;      y2 = y1 - self.thickness
      qpoints_list= [    #adding the seven points to contruct the arrow. starting from top left corner (when start mode)
                 QtCore.QPointF( x1,y4 ),  QtCore.QPointF( x3,y4 ),   QtCore.QPointF( x3,y3 ),   QtCore.QPointF( x2,y3 ),   \
                 QtCore.QPointF( x2,y2 ),  QtCore.QPointF( x3,y2 ),  QtCore.QPointF( x3,y1 ),   QtCore.QPointF( x1,y1 )     ]
      ear=   QtGui.QGraphicsPolygonItem(   QtGui.QPolygonF(  qpoints_list ),   self.item)  
      ear.setPen(  QtGui.QPen(QtCore.Qt.NoPen)  )
      ear.setBrush(fill_brush)      
      
class GeneFace(faces.Face):
  """Mother class for some Face representations created for genes. It is not functional by itself, it requires some functions (e.g. update_items) to be defined in the subclasses.
     Initialized as:
     gf= GeneFace(gene_object,   **key_attributes )
     The gene object is of class MMlib.gene; in its minimal form, it must contain a ".strand" attribute (+ or -), a .chromosome attribute, and an .id attribute. By default id and chromosome are printed within the representation, but this can be controlled with the following keyattribute: 

    *printed: a hash with as keys, some keywords of attributes of the self.gene object that will be printed as text.
      - possible values: ['id', 'chromosome', 'boundaries', 'text']  OR any other attribute found in self.gene; the vertical area will be partitioned in an equal number of horizontal stripes to accomodate all the information, with order defined by attribute order_text (any additional attribute, not included in order_text, goes as last). 'id', 'chromosome' and 'text' are looked for in the self.gene object as text attributes (note: if self.gene.text is not defined, an error occurs). For boundaries, there's an hardcoded routine that looks in the self.gene.exon coordinates and print a summary of it, using an text separator the self.gene.strand. Example: strand:+ positions:3-30,500-785   --> summary is 3+785. 
    *colors: this hash defines the colors for all objects displayed. As keys, you may provide annotation names (for subclass ExonViewFace) or any of these reserved keywords:
      - bkg:       the fill color of the main gene rectangle
      - outline:   the line color around the main gene rectangle
      - gap:       the color of the line denoting gaps  (which in dotted style)  #only for ExonViewFace
      - text:      the color of all text displayed 
"""
  default_attributes={
      'order_text':['id', 'chromosome', 'boundaries', 'text'],
      'height':40, 'width':150,  #'width_per_position':None, 
      'font':'Courier New', 'font_size':12, 'shrink_font_to_fit':True, 'font_margin':4,
      'colors':{'bkg':'#4444AA', 'outline':'#000000', 'text':'#000000' },
      'printed':{'id':1, 'chromosome':1, 'text':0, 'boundaries':0},
      'pen_size':1,
      'rotable':False,
    }

  def _width(self):        return self.item.rect().width()
  def _height(self):       return self.item.rect().height()

  def __init__(self, gene_obj, **attributes):
      faces.Face.__init__(self)
      self.type = "item"
      self.item = None
      self.gene= gene_obj
      new_attribute_keys ={}   # to keep track of those which are not defined in the default dict;  key: attribute name --> value: list of keys
      for attribute_name in attributes: 
        if self.default_attributes.has_key(attribute_name):
          if type( attributes[attribute_name] ) == dict:
            for k in attributes[attribute_name].keys():  #### just checking
              if not self.default_attributes[attribute_name].has_key(k):
                if  not (attribute_name in ['colors', 'printed']): printerr('init GeneFace WARNING  the key "'+str(k)+'" of attribute "'+attribute_name+'" is not present for this class!', 1)
                else:
                  if not new_attribute_keys.has_key(attribute_name): new_attribute_keys[attribute_name]=[]
                  new_attribute_keys[attribute_name].append(k)                         
        else:        printerr('init GeneFace WARNING  the attribute "'+attribute_name+'" is not present for this class!', 1 )
      ### parsing attributes: set
      for attribute_name in self.default_attributes:
        if attributes.has_key(attribute_name):
          if type(self.default_attributes[attribute_name]) != dict:   setattr(self, attribute_name, attributes[attribute_name])
          else: 
            h={}
            keys_to_parse=self.default_attributes[attribute_name].keys()
            if new_attribute_keys.has_key(attribute_name):  keys_to_parse += new_attribute_keys[attribute_name]
            for k in keys_to_parse:
              if  attributes[attribute_name].has_key(k):    h[k]=attributes[attribute_name][k]
              else:                                         h[k]=self.default_attributes[attribute_name][k]  
            setattr(self, attribute_name, h)            
        else:          setattr(self, attribute_name, self.default_attributes[attribute_name])  # value not defined in __init__

  def sort_fields_to_print(self, field_list):
      """ Used to determine the vertical order of things that must be printed (self.text attribute)""" 
      k_hash={}
      for field in field_list:
        position=None
        #finding index of field in self.order_text
        for index, field_name in enumerate(self.order_text): 
          if field==field_name: 
            position=index
            break
        if position is None: position=len(self.order_text)
        k_hash[field]=position
      return sorted(field_list, key=lambda x:k_hash[x])

  def update_items(self):        raise Exception, "ERROR calling the raw update_items function from class GeneView! This function must be defined in the subclasses. Aborting! "

  def draw_text_fields(self):
      ### printing any text 
      fields_to_print=self.sort_fields_to_print( [ f for f in self.printed if self.printed[f] ] )  # margin #--area-field1--# margin #--area_field2--# margin ... #--area_fieldN--# margin #
      if fields_to_print:
        area_field_height=( self.height - self.font_margin*(len(fields_to_print)+1) ) / float( len(fields_to_print) )
        font_used=  QtGui.QFont( self.font, self.font_size )    
        brush_used= QtGui.QColor( self.colors['text']   )
        for field_index, field_name in enumerate( fields_to_print ):
          text_obj = QtGui.QGraphicsSimpleTextItem()
          text_obj.setBrush( brush_used  )
          text_obj.setFont( font_used )
          if field_name.startswith('&'):      the_text = str(  getattr(self.gene, field_name[1:]) () )  #it's a function; written like &function_name
          elif field_name == 'boundaries':    
            join_symbol=self.gene.strand
            the_text = join(map(str, self.gene.boundaries()), join_symbol)
          else:                     the_text = str(  getattr(self.gene, field_name)  )  #generic for attribute
          text_obj.setText(the_text)
          if self.shrink_font_to_fit:  
            shrink_font_to_box(text_obj, self.width, area_field_height, category=field_name)  #reducing font if necessary to fit the allocated space
      ## move position to center it
          br=text_obj.boundingRect();       x_text= ( self.width - br.width()) /2;    
          y_text= self.font_margin +   (area_field_height+self.font_margin)*field_index   +   ( area_field_height - br.height() )/2    
          text_obj.setPos( x_text, y_text )
          text_obj.setParentItem(self.item)

class GeneViewFace(GeneFace):
  """ Typical GFF visualization, representing exons and introns. Something like this:
    _______     ________
 --|       |---|        |--     
    -------     --------
Usage:     gvf= GeneViewFace(gene_object, introns=[], xmin=None, xmax=None, shrink_introns=None, **key_attributes )
     The gene object is of class MMlib.gene; in its minimal form, it must contain a ".strand" attribute (+ or -), a .chromosome attribute, and an .id attribute. 
    *introns can be used to display connectors in the typical triangle style. These can link positions which are not necessarily exon boundaries.  Argument must be provided in the form: introns=[ [start1,stop1], [start2, stop2] .. ], where start < stop
    *xmin and *xmax define the boundaries of the genomic region displayed. If None, the min and max values present in the gene exons or in the introns argument are used
    *shrink_introns: often introns are so much larger than exons that it is bad for visualization. Set this to any integer number to set the maximum intron size which is faithfully represented; introns bigger than this will be displayed as big as this, with a pattern indicating the compression.
     By default the .id attribute of the gene object (and only that) is printed within the representation, but this can be controlled with the following keyattribute: 
    *printed: a hash with as keys, some keywords of attributes of the self.gene object that will be printed as text.
      - possible values: ['id', 'chromosome', 'boundaries', 'text']  OR any other attribute found in self.gene; the vertical area will be partitioned in an equal number of horizontal stripes to accomodate all the information, with order defined by attribute order_text (any additional attribute, not included in order_text, goes as last). 'id', 'chromosome' and 'text' are looked for in the self.gene object as text attributes (note: if self.gene.text is not defined, an error occurs). For boundaries, there's an hardcoded routine that looks in the self.gene.exon coordinates and print a summary of it, using an text separator the self.gene.strand. Example: strand:+ positions:3-30,500-785   --> summary is 3+785. 
    *colors: this hash defines the colors for all objects displayed. As keys, you may provide annotation names (for subclass ExonViewFace) or any of these reserved keywords:
      - bkg:       the fill color of the main gene rectangle
      - outline:   the line color around the main gene rectangle
      - text:      the color of all text displayed 
      - line:      the color of the genome line
      - intron:    the color of the introns drawn
      - exonN or exonN.outline:   define a key like this, with N replacing the index of any exon (1 is first at 5'), to change the color of the background (or outline) of any specific exon. Useful to highlight an exon
    SIZE:
     for height, you must define either: *height (total height of face in pixels) or *height_exons and height_introns (to define manually the space attributed to the exon and the intron linker.
     Similarly for width you must define either *width (total width of face in pixels)  or *width_per_position, where position refers to the coordinates in the gene object
    *pen_size defines the size of the pen used to draw all lines (default: 1)
    *intron_pen_size defines the size of the pen used for introns only (default: 2)
"""
  default_attributes={
      'order_text':['id', 'chromosome', 'boundaries', 'text'],
      'height':None, 'height_exons':40,  'height_introns':15, 
      'width':150,  'width_per_position':None, 
      'font':'Courier New', 'font_size':12, 'shrink_font_to_fit':True, 'font_margin':4,
      'colors':{'bkg':'#4444AA', 'outline':'#000000', 'text':'#000000', 'line':'#000000', 'intron':'#000000', #'exon1':'#888844', 'exon1.outline':'#AA0000'
      },
      'printed':{'id':1, 'chromosome':0, 'text':0, 'boundaries':0},
      'pen_size':1,       'intron_pen_size':2,
      'rotable':False,
    }

  def __init__(self, gene_obj, introns=[], xmin=None, xmax=None,  shrink_introns=None, **attributes):
      GeneFace.__init__(self, gene_obj, **attributes)      ### attributes are added / updated here
      ## setting height
      if not introns: self.height_introns=0
      if self.height is None:    self.height = self.height_exons  +   self.height_introns  
      elif not introns:          self.height_exons = self.height
      else:             
        self.height_exons =    int(0.6 * self.height)
        self.height_introns = self.height - self.height_exons

      # to determine self.width_per_position I need to know how many introns are compressed. 
      self.compressed_regions=[]   #[  [start, end, shrinked_to], ...  ]    --> shrinked_to meaning, the apparent width in the compressed version
      self.compressed_intron_indexes={}
      compressed_to_length    = 0;         total_length_compressed = 0
      if shrink_introns:
        for intron_index, intron_positions in enumerate( self.gene.introns().exons ):
          intron_start, intron_end = intron_positions
          if intron_end - intron_start +1 > shrink_introns:
            self.compressed_regions.append(  [intron_start, intron_end, shrink_introns]  )
            self.compressed_intron_indexes[intron_index]=True
        total_length_compressed = sum(map( lambda x:x[1]-x[0]+1, self.compressed_regions))
        compressed_to_length  = shrink_introns * len(self.compressed_regions)
                                    
      ### get or set self.width and self.width_per_position   and min and max boundaries
      boundaries=self.gene.boundaries()
      if xmin is None: xmin=min( boundaries +  [i_start for i_start,i_end in introns] + [i_end for i_start,i_end in introns] )
      if xmax is None: xmax=max( boundaries +  [i_start for i_start,i_end in introns] + [i_end for i_start,i_end in introns] )
      if self.width_per_position:     self.width=             self.width_per_position *  (xmax-xmin+1 - total_length_compressed + compressed_to_length )
      elif self.width:                self.width_per_position=float(self.width) /        (xmax-xmin+1 - total_length_compressed + compressed_to_length )
      else: raise Exception, "GeneViewFace ERROR neither width or width_per_position are set!"

      self.processed_coords = {'exons':[], 'introns':[]} #, 'compressed':[]}   # [ [x1, x2] ...  ]  in pixels   ##### referring to exon positions as "positions" and to pixel positions as "coords"
      for exon_index, exon_positions in enumerate(self.gene.exons):
        exon_start, exon_end= exon_positions
        if self.gene.strand=='-':                     delta_pos_start = xmax - exon_end;          delta_pos_end   = xmax - exon_start          
        else:                                         delta_pos_start = exon_start - xmin;        delta_pos_end   = exon_end   - xmin 

        ## correcting for compressed introns
        for compressed_index, compr_indexes in enumerate( self.compressed_regions ):
            compressed_start, compressed_end, shrinked_to = compr_indexes
            all_compressed_are_behind= exon_start < compressed_start
            if self.gene.strand=='-': all_compressed_are_behind = not all_compressed_are_behind
            if all_compressed_are_behind:  break
            diff_for_compression = shrinked_to - (compressed_end - compressed_start +1) 
            delta_pos_start += diff_for_compression;  delta_pos_end += diff_for_compression
                    
        coord_start = delta_pos_start * self.width_per_position;   coord_end = delta_pos_end * self.width_per_position 
        self.processed_coords['exons'].append([coord_start, coord_end])

      for intron_index, intron_positions in enumerate(introns):
        intron_start, intron_end= intron_positions
        if self.gene.strand=='-':          delta_pos_start = xmax - intron_end;         delta_pos_end   = xmax - intron_start  
        else:                              delta_pos_start = intron_start - xmin;       delta_pos_end   = intron_end   - xmin 
        ## correcting for compressed introns
        #print "******", self.gene.id,  intron_start, intron_end
        for compressed_index, compr_indexes in enumerate( self.compressed_regions ):
            compressed_start, compressed_end, shrinked_to = compr_indexes
          #  print compressed_start, compressed_end, shrinked_to
            if    self.gene.strand=='+':  all_compressed_are_behind=  intron_end   < compressed_start  # if this intron is also compressed, this becomes False, but we add it to end only 
            elif  self.gene.strand=='-':  all_compressed_are_behind=  intron_start > compressed_end

            diff_for_compression = shrinked_to - (compressed_end - compressed_start +1) 
            if all_compressed_are_behind:                break
            #else:               print "summing diff!", compressed_index, diff_for_compression
            if    (self.gene.strand=='+' and intron_start>compressed_end) or (self.gene.strand=='-' and intron_end<compressed_start):
              delta_pos_start += diff_for_compression
         #     print "diff start", diff_for_compression
         #   print "diff end", diff_for_compression
            delta_pos_end += diff_for_compression

        coord_start     = delta_pos_start * self.width_per_position;   coord_end     =  delta_pos_end * self.width_per_position
        self.processed_coords['introns'].append([coord_start, coord_end])
      self.xmin=xmin;       self.xmax=xmax; self.introns=introns;     self.shrink_introns=shrink_introns
                
  def update_items(self):
      """ the actual function doing the graphics """
      self.item = QtGui.QGraphicsRectItem(0, 0, self.width, self.height) # backbone object, no color 
      self.item.setPen(QtGui.QPen(QtCore.Qt.NoPen)) #no black border      

      if self.colors['outline'] is None: outline_color= self.colors['bkg']
      else: outline_color= self.colors['outline']
      outline_pen = QtGui.QPen(    QtGui.QColor(outline_color),  self.pen_size, QtCore.Qt.SolidLine)
      exon_brush  = QtGui.QBrush(  QtGui.QColor(self.colors['bkg']) )
      line_pen    = QtGui.QPen(    QtGui.QColor( self.colors['line']  ),  self.pen_size, QtCore.Qt.SolidLine)
      intron_pen =  QtGui.QPen(    QtGui.QColor(self.colors['intron']),   self.intron_pen_size, QtCore.Qt.SolidLine)

      if self.processed_coords['exons']:
        ## draw line up to first exon box
        coord_start_line = 0
        coord_end_line   = self.processed_coords['exons'][0][0]-1
        if int( coord_end_line - coord_start_line  ) > 0:
#          print 'drawing first line', coord_start_line, coord_end_line
          line =QtGui.QGraphicsLineItem( coord_start_line   , self.height_introns + self.height_exons / 2.0,   coord_end_line,   self.height_introns + self.height_exons / 2.0,           self.item)             
          line.setPen(line_pen);          
        ####### drawing exons
        for exon_index, exon_coords in enumerate(self.processed_coords['exons']):
          coord_start, coord_end = exon_coords
          colored_rect = QtGui.QGraphicsRectItem(  coord_start, self.height_introns, coord_end-coord_start+1, self.height_exons, self.item)
          ### custom color specified for this exon?
          if  'exon'+str(exon_index+1) in self.colors:   colored_rect.setBrush(QtGui.QBrush(QtGui.QColor(self.colors['exon'+str(exon_index+1)])))
          else:           colored_rect.setBrush(exon_brush)          
          if  'exon'+str(exon_index+1)+'.outline' in self.colors:   colored_rect.setPen(QtGui.Pen(QtGui.QColor(self.colors['exon'+str(exon_index+1)+'.outline']), self.pen_size, QtCore.Qt.SolidLine  ))
          else:           colored_rect.setPen(outline_pen)

          ## drawing line just after this
          if   exon_index+1 < len(  self.processed_coords['exons']  ):  ## not last exon
            coord_start_line=coord_end +1 
            coord_end_line=  self.processed_coords['exons'][exon_index+1][0] -1   #start of next exon
            #print "drawing ", coord_end+1, coord_next_exon_start-1
          else:  # last exon
            coord_start_line=coord_end +1 
            coord_end_line  = self.width
          if int( coord_end_line - coord_start_line  ) > 0:
            #print "drawing last line", coord_start_line, coord_end_line
            if not exon_index in self.compressed_intron_indexes:
              line =QtGui.QGraphicsLineItem( coord_start_line   , self.height_introns + self.height_exons / 2.0,   coord_end_line,   self.height_introns + self.height_exons / 2.0,      self.item);                 line.setPen(line_pen);          
            else:
              compressed_start, compressed_end = coord_start_line, coord_end_line
              size= compressed_end - compressed_start +1

        ####### drawing compressed      
        ##          p4 / /p7              y2
        ##   p1     p2/ /p6
        ##    -------/ /------- p8        y1
        ##          / /
        ##       p3/ /p5                  y0
        ##
              y0=self.height_introns + 0.9*self.height_exons 
              y1=self.height_introns + self.height_exons/2.0
              y2=self.height_introns + 0.1*self.height_exons 
              p1=compressed_start
              p2=compressed_start + 0.4 * size
              p3=compressed_start + 0.3 * size
              p4 = p5 = compressed_start + 0.5 * size
              p6=compressed_start + 0.6 * size
              p7=compressed_start + 0.7 * size
              p8= compressed_end          
              # horiz line #1
              line =QtGui.QGraphicsLineItem( p1, y1, p2, y1, self.item);         line.setPen(line_pen);
              # horiz line #2
              line =QtGui.QGraphicsLineItem( p6, y1, p8, y1, self.item);         line.setPen(line_pen);
              # rising line #1
              line =QtGui.QGraphicsLineItem( p3, y0, p4, y2, self.item);         line.setPen(line_pen);
              # rising line #2
              line =QtGui.QGraphicsLineItem( p5, y0, p7, y2, self.item);         line.setPen(line_pen);

      if self.processed_coords['introns']:        
        ####### drawing introns
        for intron_index, intron_coords in enumerate(self.processed_coords['introns']):
          coord_start, coord_end = intron_coords
          halfway= (coord_end + coord_start) / 2.0
          # rising line
          line =QtGui.QGraphicsLineItem( coord_start, self.height_introns , halfway, 0 , self.item)
          line.setPen(intron_pen);
          # descending line
          line =QtGui.QGraphicsLineItem( halfway, 0 , coord_end, self.height_introns , self.item)
          line.setPen(intron_pen);
        
      self.draw_text_fields()
      return  

class ArrowGeneFace(GeneFace):
  """ Thought for syntheny; each gene is displayed as a colored arrow, towards right or left (+ or - strand). If the gene.strand is None, it's just a rectangle (no arrow)
Initialize like:     agf= ArrowGeneFace( gene_object,  **key_attributes)
Gene object is of class MMlib.gene; these attributes must be defined: strand, chromosome, id.
Example to obtain the gene object:   gene_object=gene(strand='+', chromosome='chr1', id='Flybase0001')
Possible **key_attributes:
  -height, width: these refer to the (invisible) external box surrounding the arrow, in pixels
  -arrow_height and pointer_width ; see figure below
  -colors: a hash defining the colors to draw this. Example: colors={'bkg':'lightblue', 'outline':'#000000', 'text': 'red'}.  See GeneFace __doc__ string for more details. If outline is None, the same color as for bkg is used. If box_bkg or box_line are not provided, none (transparent) is used
  -printed: decides what text is printed in the figure. Text is provided as attributes of the gene object. See GeneFace __doc__ string for more details. By default in this class, only the gene.id attribute is printed
  -others: rotable, font, font_size, shrink_font_to_fit, font_margin : see GeneFace

     < ---      width      --- >
  ^  ###########################
  |  #  <box_bkg>        |\    #
  |  #-------------------|  \  #   -
 hei #|      < bkg >         ) #   |--> arrow_height
 ght #-------------------|  /  #   -
  |  #      outline ^    |/    #
  -  ###########################
     box_line ^          <--->
                      pointer_width
  """
  default_attributes={
      'order_text':['id', 'chromosome', 'text', 'boundaries'],
      'height':40, 'width':150,  #'width_per_position':None, 
      'arrow_height':25,      'pointer_width':30,
      'font':'Courier New', 'font_size':8, 'shrink_font_to_fit':False, 'font_margin':4,
      'colors':{'bkg':'#4444AA', 'outline':None, 'text':'#000000', 'box_bkg':None, 'box_line':None },
      'printed':{'id':1, 'chromosome':0, 'text':0, 'boundaries':0},
      'pen_size':1,
      'rotable':False,
    }
  
  def update_items(self):
      """ the actual function doing the graphics """
      if not self.gene.strand in ['+', '-', None]: raise Exception, "ArrowGeneFace ERROR this face must be initialized with a gene object with a .strand attribute which is + or -  (or None) ! This was received: {0}".format(self.gene.strand)
      self.item = QtGui.QGraphicsRectItem(0, 0, self.width, self.height) # backbone object, no color 
      if not 'box_line' in self.colors or self.colors['box_line'] is None:   box_pen=QtGui.QPen(QtCore.Qt.NoPen)
      else:                                                                  box_pen=QtGui.QPen(  QtGui.QColor(self.colors['box_line']), self.pen_size, QtCore.Qt.SolidLine )
      if not 'box_bkg' in self.colors or self.colors['box_bkg'] is None:     box_brush=QtGui.QBrush(  QtCore.Qt.NoBrush  )
      else:                                                                  box_brush=QtGui.QBrush(  QtGui.QColor(self.colors['box_bkg'])  )
      self.item.setPen(box_pen) #no black border      
      self.item.setBrush(box_brush)      

      if not 'outline' in self.colors or self.colors['outline'] is None: outline_color= self.colors['bkg']
      else:           outline_color= self.colors['outline']
      outline_pen = QtGui.QPen(    QtGui.QColor(outline_color),  self.pen_size, QtCore.Qt.SolidLine)
      rect_brush  = QtGui.QBrush(  QtGui.QColor(self.colors['bkg']) )
      #rect_pen =   QtGui.QPen(    QtGui.QColor(self.colors['bkg']),  self.pen_size, QtCore.Qt.SolidLine)

      y1_arrow_trunk=  (self.height - self.arrow_height) / 2.0
      y2_arrow_trunk=  y1_arrow_trunk + self.arrow_height
      y1_arrow_head =  0
      y2_arrow_head =  self.height
      y_arrow_head=    self.height/2.0
      arrow_trunk_width= self.width - self.pointer_width
      x1_arrow_head=  arrow_trunk_width
      x2_arrow_head=  self.width
      x_arrow_trunk = 0
        
      if self.gene.strand=='-':
        y1_arrow_trunk, y2_arrow_trunk= y2_arrow_trunk, y1_arrow_trunk
        y1_arrow_head, y2_arrow_head  = y2_arrow_head, y1_arrow_head
        x_arrow_trunk = self.width
        x1_arrow_head = self.pointer_width
        x2_arrow_head = 0
      
      qpoints_list= [    #adding the seven points to contruct the arrow. starting from top left corner (when plus strand)
                 QtCore.QPointF(   x_arrow_trunk,  y1_arrow_trunk)  , #A
                 QtCore.QPointF(   x_arrow_trunk,  y2_arrow_trunk)  , #B
                 QtCore.QPointF(   x1_arrow_head,  y2_arrow_trunk)  , #C
                 QtCore.QPointF(   x1_arrow_head,  y2_arrow_head)   , #D
                 QtCore.QPointF(   x2_arrow_head,  y_arrow_head)    , #E -- pointer
                 QtCore.QPointF(   x1_arrow_head,  y1_arrow_head)   , #F
                 QtCore.QPointF(   x1_arrow_head,  y1_arrow_trunk)  ] #G

      if self.gene.strand is None:   ## just a rectangle
        qpoints_list = [ QtCore.QPointF(x_arrow_trunk,  y1_arrow_trunk), QtCore.QPointF(x_arrow_trunk,  y2_arrow_trunk),  
                         QtCore.QPointF(x2_arrow_head,  y2_arrow_trunk), QtCore.QPointF(x2_arrow_head,  y1_arrow_trunk)   ]

      arrow=   QtGui.QGraphicsPolygonItem(   QtGui.QPolygonF(  qpoints_list ),   self.item)  
      arrow.setPen(outline_pen)
      arrow.setBrush(rect_brush)      
      self.draw_text_fields()
      return  
 
def syntheny_view ( list_of_genes, get_color= lambda x:x.color, face=ArrowGeneFace, margins={'margin_left':4, 'margin_right':4, 'margin_top':4, 'margin_bottom':4 }, **other_face_attributes):
  """  Function that accepts a list of genes in input (MMlib.gene class) indicating a synthenic block  and returns a list of ArrowGeneFace faces to represent this.
To represent things other than genes, you can include any of these special values as elements in the input list of genes:
   '-'      :  gaps in the syntheny  (-> GapFace)
   ''       :  white spacer          (-> GapFace  with color=None)
   'start'  :  contig starts here    (-> ContigTerminusFace with mode='start')
   'start'  :  contig starts here    (-> ContigTerminusFace with mode='start')
This method should be called for each node in which you want to have this view. The representation of genes in different node can be conceptually linked through colors.
The function get_color is a generalized way to derive the color of each gene. This function should return a color name.
Margins is a hash defining the margin attributes set in the faces created.
  other_face_attributes:  other keywords can be specified setting any attribute of the face (e.g. see ArrowGeneFace.default_attributes). All attributes specified like this are set for all genes in list_of_genes.
  NOTE: it is also possible to specify attributes in the input gene objects, allowing to customize them one by one. The attribute "face_attributes" must be defined in the gene object; this is a hash in the same format as you would put other_face_attributes. 
   E.g. g=gene(strand='+', chromosome='chr1', face_attributes ={'font':'Arial', 'font_size':5}  )    --> then provide g in list_of_genes
Additionally, four special attributes in the gene objects are recognized: color_outline, color_text, color_box_bkg, color_box_line ; if these are defined, the default values are overriden.
       
Usage example: 
for col_index, the_face in enumerate( syntheny_view( [ gene(strand='+', id='gene1',   color='blue', color_box_bkg='yellow'), gene(strand='-', id='gene2', color='green',     color_outline='black'  )])):   node.add_face( the_face, col_index+1, 'aligned')

 """
  #face_container= FaceContainer()
  face_list=[]
  for g_index, g in enumerate(list_of_genes):
    if 'width' in  other_face_attributes: width=other_face_attributes['width']
    else:                                 width=face.default_attributes['width'] 
    if 'height' in  other_face_attributes: height=other_face_attributes['height']
    else:                                  height=face.default_attributes['height'] 

    if      g is None or g=='-':  agf=GapFace(width=width, height=height)
    elif    g == '':              agf=GapFace(width=width, height=height, color=None)
    elif g in ['start', 'end']:   agf=ContigTerminusFace(width=width, height=height, mode=g) 
    else:
      try:     
        bkg_color= get_color( g )
        assert bkg_color
      except:  
        printerr("ERROR syntheny_view the function get_color provided failed with this gene: "+str(g.id), 1)
        raise

      if hasattr(g, 'color_outline') and g.color_outline: outline_color=g.color_outline
      else:                                               outline_color=bkg_color
      if hasattr(g, 'color_text') and g.color_text: text_color=g.color_text
      else:                                         text_color='black'
      if hasattr(g, 'color_box_bkg')  and g.color_box_bkg:   box_bkg_color=g.color_box_bkg
      else:                                                  box_bkg_color= None
      if hasattr(g, 'color_box_line') and g.color_box_line:  box_line_color=g.color_box_line
      else:                                                  box_line_color=None

      colors={'bkg':bkg_color, 'outline':outline_color, 'text':text_color, 'box_bkg':box_bkg_color, 'box_line':box_line_color  }      
      #print colors
      
      other_face_attributes_for_this= dict(other_face_attributes) #copying dict
      if hasattr(g, 'face_attributes'): 
        for k in g.face_attributes:          
          if    k=='width':  width= g.face_attributes[k]
          elif  k=='height': height=g.face_attributes[k]
          elif  k=='colors': colors=g.face_attributes[k]
          else: other_face_attributes_for_this[k] = g.face_attributes[k]
      for k in ['width', 'height', 'colors']: 
        if k in other_face_attributes_for_this: del other_face_attributes_for_this[k]

      agf=face(g, colors=colors, width=width, height=height, **other_face_attributes_for_this)
    for k_margin in margins:      setattr(agf, k_margin, margins[k_margin])
    face_list.append(agf)
  return face_list


class ExonViewFace(GeneFace):
    """ Exon centric representation of a gene with potentially multiple exons, as a colored rectangle. Introns are compressed to vertical lines, and gaps are permitted so that multiple instances of this could represent a protein alignment. Attributes may be printed within the rectangle.
  -Init: ExonViewFace(   gene_obj, gaps="", printed={}, annotations={}, **attributes   )
    *gene object: must have a strand, a chromosome, and list of exon positions ( .exons attribute of gene, like [ [st, end], [st2, end2], ...  ]  )
    *gaps: three types of arguments can be provided here; 
      - a list like [ (pos_of_gap_in_nucleotide, gap_length) ]
      - a string that, removing gaps (-), has the same length as the gene object (~ aligned nucleotides); gaps are indicated with "-", the other chars are ignored
      - similarly, a string with sequence of length like gene_length/3, (~aligned proteins)
    *printed:   see GeneFace help
    *annotations:   a hash describing the information of any positional feature that you want to be displayed. format: { 'annotation_name':  {'map': 'MAP_TYPE',  'start':[3, 40, 80], 'end':[8, 46, 81]} } ; additionally: it can also have y1, and it can also have y0.   y0-y1 define the vertical space that is span by this feature. In practive it can be used to display another layer of information for each feature, like a score. Note: y1 and y0 are floats between 0 and 1. 0 means the inferior baseline, 1 the superior baseline.
      - MAP_TYPE defines the relation between the positions given and the actual drawing position; note that positions are always 1-based, nucleotide based; value "SEQ" --> the positions are relative to the single sequence (e.g. 1 is first position of gene object); value "ALI" --> the positions are relative to the alignment, meaning, gaps are also counted as possible positions (e.g. 1 is always the first position in the drawing, the extreme left of the rectangle)
      - 'start' and 'end' provide the list of positions. They should have the same number of elements. In case the annotation refers to a single point, 'end' can be omitted, and in this case it's like end= start for each of the annotation instances.
      (color of the annotations can be defined by providing through the "colors" argument of __init__: just use as key the annotation name. If annotations spans more than a single position ('end' defined) the color defines the pen outline and also the brush fill; otherwise, only pen is defined)
      the annotation names 'introns' is reserved for drawing the introns as vertical lines. The value must be True to display.   ###and 'frameshifts' and stop
    *colors:    see GeneFace help
    *other attributes: see GeneFace help; list: order_text, height, width OR width_per_position, font, font_size, shrink_font_to_fit, font_margin, rotable
"""
    default_attributes={
      'order_text':['id', 'chromosome', 'boundaries', 'text'],
      'height':40, 'width':150, 'width_per_position':None, 'min_for_gaps':5,
      'font':'Courier New', 'font_size':12, 'shrink_font_to_fit':True, 'font_margin':4,
      'colors':{'bkg':'#4444AA', 'outline':'#000000', 'text':'#000000', 'gap':'#999999', 'introns':'#000000', 'frameshifts':'#EE2222' },
      'printed':{'id':1, 'chromosome':1, 'text':0, 'boundaries':1},
      'annotations':{'introns':True, 'frameshifts':True},
      'pen_size':1,
      'rotable':False,
    }

    def __init__(self, gene_obj, gaps="", annotations=None, **attributes):
      GeneFace.__init__(self, gene_obj, **attributes)      ### attributes are added / updated here

      ###  processing gaps position, by joining them
      self.gaps=[]    #form:    [gap_relative_pos, gap_length]   relative pos in in nucleotides, 0 means before the first, 1 after the first etc          
      if type(gaps)==str and len(nogap(gaps))*3  == len(self.gene) : # if variable gaps is provided as aligned proteins
        gaps=join([ 3*char   for char in gaps ],'')  ## modifying gaps to enter in next "if"
      if type(gaps)==str and len(nogap(gaps)) == len(self.gene): # if variable gaps is provided as aligned nucleotides
        index=0; tot_gaps=0 #increasing
        while index<len(gaps):
          if gaps[index]=='-': 
            gap_start=index
            while index<len(gaps) and gaps[index]=='-':               index+=1
            gap_length=index-gap_start
            self.gaps.append([gap_start-tot_gaps, gap_length])
            tot_gaps+=gap_length            
          index+=1
      elif type(gaps)==list: self.gaps=gaps ## this was provided directly in the form it is stored (See above)
      else:         raise Exception, "ERROR initializing ExonViewFace (id="+self.gene.id+") the gaps argument was not recognized! wrong length maybe? \ngaps="+str(gaps)+'\nlen(gaps)='+str(len(gaps))+' gene.length()='+str(self.gene.length())
      # add gaps of length 0 at beginning and end, if not present; this allows to have a nice data structure for later parsing
      if not self.gaps or self.gaps[0][0] != 0:   self.gaps.insert(0, [0, 0])  
      if self.gaps[-1][0] != self.gene.length():  self.gaps.append([self.gene.length(), 0])        

      ### get or set self.width and self.width_per_position
      if self.width_per_position:     self.width= self.width_per_position * (  self.gene.length() + sum([gap_length for g_p, gap_length in self.gaps])  )
      elif self.width:                self.width_per_position=float(self.width)/ ( self.gene.length() + sum([gap_length for g_p, gap_length in self.gaps]) )
      else: raise Exception, "ExonViewFace ERROR neither width or width_per_position are set!"

      ### filtered gaps: precomputed minimal data structure on gaps positiions and length
      self.filtered_gaps=[]  #pos, length, gap_length_before_this   #### I need to keep memory of the total size of all gaps -- including those too short to be displayed
      total_gap_length=0   ## in nucleotides, progressively increased      
      min_length_gap_to_be_displayed= self.min_for_gaps / self.width_per_position
      for pos, length in self.gaps:
        if (length >=min_length_gap_to_be_displayed) or (pos in [0, self.gene.length() ]) :
          self.filtered_gaps.append(  [pos, length, total_gap_length]  )
        total_gap_length+=length

      #### PARSE ANNOTATIONS; make everything nucleotide, alignment based
        # preparing "annotations";  filling default values 
      if annotations is None: annotations={}
      for key in  self.default_attributes['annotations']:  
        if not annotations.has_key(key): annotations[key]= self.default_attributes['annotations'][key]
      #  print "pre: ", annotations, '\n', self.annotations
      ## pre processing special annotation: introns
      if (annotations.has_key('introns') and annotations['introns']) or (annotations.has_key('frameshifts') and annotations['frameshifts']):
        introns_dict=    {'map':'SEQ', 'start':[], 'end':[]}
        frameshifts_dict={'map':'SEQ', 'start':[], 'end':[]}
        if len(self.gene.exons)>=2:  
          nucleotides_behind=0
          for exon_index in range( 1, len(self.gene.exons) ):
            previous_exon_length= self.gene.exons[exon_index-1][1] - self.gene.exons[exon_index-1][0]
            if    self.gene.strand=='+':    intron_length= self.gene.exons[exon_index][0] - self.gene.exons[exon_index-1][1] - 1 
            elif  self.gene.strand=='-':    intron_length= self.gene.exons[exon_index-1][0] - self.gene.exons[exon_index][1] - 1 
            nucleotides_behind+=previous_exon_length
            if intron_length < 6 : 
              frameshifts_dict['start'].append( nucleotides_behind )
              frameshifts_dict['end'].append( nucleotides_behind )
            else: 
              introns_dict['start'].append( nucleotides_behind )
              introns_dict['end'].append( nucleotides_behind )
        if annotations.has_key('introns') and annotations['introns']  and introns_dict['start']:   #it was requested, and there's at least one intron to display
          annotations['introns']= introns_dict
        elif   annotations.has_key('introns'):           del annotations['introns']          

        if annotations.has_key('frameshifts') and annotations['frameshifts']  and frameshifts_dict['start']:   #it was requested, and there's at least one fs to display
          annotations['frameshifts']= frameshifts_dict
        elif annotations.has_key('frameshifts'):           del annotations['frameshifts']          
                       
      ## general parse annotations; putting everything in alignment relative coordinates; generating and filling self.annotations
      self.annotations={}
      for annotation_name in annotations:
        if not annotations[annotation_name]: 
          self.annotations[annotation_name]=False
          continue
        ann_dict = annotations[annotation_name]
        if ann_dict.has_key('end') and len(ann_dict['start'])!= len(ann_dict['end']): raise Exception, "ExonViewFace ERROR annotation name: "+annotation_name+" start and end index lists have different length! start: "+str(ann_dict['start'])+' end: '+str(ann_dict['end'])
        elif not ann_dict.has_key('end'):  ann_dict['end']=ann_dict['start']
        ### checking y1 and y0
        if not ann_dict.has_key('y1'): ann_dict['y1']=  [1.0 for index in range( len(ann_dict['start']) ) ]
        elif len(ann_dict['start'])!= len(ann_dict['y1']): raise Exception, "ExonViewFace ERROR annotation name: "+annotation_name+" start and y1 values lists have different length! start: "+str(ann_dict['start'])+' y1: '+str(ann_dict['y1'])
        if not ann_dict.has_key('y0'): ann_dict['y0']=  [0.0 for index in range( len(ann_dict['start']) ) ]
        elif len(ann_dict['start'])!= len(ann_dict['y0']): raise Exception, "ExonViewFace ERROR annotation name: "+annotation_name+" start and y0 values lists have different length! start: "+str(ann_dict['start'])+' y0: '+str(ann_dict['y0'])          

        if any( [  y0_value<0.0 or y0_value>1.0  for y0_value in ann_dict['y0'] ] ): 
          raise Exception, "ExonViewFace ERROR annotation name: "+annotation_name+" y0 values not valid! they must be between 0.0 and 1.0! "+str(ann_dict['y0'])
        if any( [  y1_value<0.0 or y1_value>1.0  for y1_value in ann_dict['y1'] ] ):
          raise Exception, "ExonViewFace ERROR annotation name: "+annotation_name+" y1 values not valid! they must be between 0.0 and 1.0! "+str(ann_dict['y1'])
        
        ### here creating self.annotations[annotation_name] !!!!
        if ann_dict['map']=='ALI':            self.annotations[annotation_name]={'start':ann_dict['start'], 'end':ann_dict['end'],  'y0':ann_dict['y0'], 'y1':ann_dict['y1']  }
        elif ann_dict['map']=='SEQ':          self.annotations[annotation_name]={'start': map(self.position_in_ali, ann_dict['start']), 'end':map(self.position_in_ali, ann_dict['end']), 'y0':ann_dict['y0'], 'y1':ann_dict['y1'] }
        else: raise Exception, "ExonViewFace ERROR map argument not recognized: '"+ann_dict['map']+"'"
        ## checking annotation positions
        ali_tot_length=self.ali_length()
        for index, start in enumerate(self.annotations[annotation_name]['start']):
          end=self.annotations[annotation_name]['end'][index]          
          if max( (start, end) ) > ali_tot_length:
            raise Exception, "ExonViewFace ERROR annotation_name: "+annotation_name+" annotation_index: "+str(index)+ " Feature out of boundaries! Total positions: "+str(ali_tot_length)+'  Start, end for this annotation were: '+str((start, end))
        
      #print   "# !! ", self.gene.id, "self.annotations", self.annotations
          
    def position_in_ali(self, pos_in_seq):
      """ Analogous to the function in class alignment; parses self.gaps to translate a sequence relative position to a alignment postion. Evrything is one based"""
      if type(pos_in_seq)!=int: raise Exception, "ERROR id "+self.gene.id+' position_in_ali called with something else than a integer: '+str(pos_in_seq)+' ('+str(type(pos_in_seq))+')'
      gaps_behind=0; gap_index=0; 
      #print [pos_in_seq], self.gaps
      while  self.gaps and  self.gaps[gap_index][0]  < pos_in_seq:
          #print "ok: ", pos_in_seq, gap_index, self.gaps[gap_index], self.gaps[gap_index][0], " < ", pos_in_seq, self.gaps[gap_index][0]  < pos_in_seq
        gaps_behind+=self.gaps[gap_index][1]
        gap_index+=1
      return pos_in_seq+gaps_behind

    def ali_length(self):
      return self.gene.length() + sum(  [g_length   for g_p, g_length in self.gaps]  )

        
    def update_items(self):
      """ the actual function doing the graphics """
      self.item = QtGui.QGraphicsRectItem(0, 0, self.width, self.height) # backbone object, no color 
      self.item.setPen(QtGui.QPen(QtCore.Qt.NoPen)) #no black border      

      gap_pen     = QtGui.QPen(    QtGui.QColor(self.colors['gap']),      self.pen_size, QtCore.Qt.DotLine)
      outline_pen = QtGui.QPen(    QtGui.QColor(self.colors['outline']),  self.pen_size, QtCore.Qt.SolidLine)
      rect_brush  = QtGui.QBrush(  QtGui.QColor(self.colors['bkg']) )
      rect_pen    = QtGui.QPen(    QtGui.QColor(self.colors['bkg']),      self.pen_size, QtCore.Qt.SolidLine)           

      for gap_index, gap_item in enumerate(self.filtered_gaps): #between a gap and the next one   ### artificial gaps are present here at pos 0 and pos length.
        gap_relative_pos, gap_length, total_gap_length= gap_item  ## relative pos is the pos of the gap in the target, nucleotide based; 0 means before the first nucleotide. 1 after the first nucleotide, etc             
        ### draw gap
        if gap_length: ## in filtered_gaps, we have either gaps of an accepted (>threshold) length, or gaps of length 0, as first or last element, just for structure
          x_start=(total_gap_length + gap_relative_pos)*self.width_per_position
          x_end=   x_start + self.width_per_position*gap_length -1        # -1 to avoid overlaps                     
          if x_end > x_start:   # basically all the times. But if we have a very small gap_length at the beginning, this could give a negative x_end, so we have this control
            line_above=QtGui.QGraphicsLineItem( x_start, 0,           x_end, 0,           self.item) 
#          line_below=QtGui.QGraphicsLineItem( x_start, self.height, x_end, self.height, self.item) 
            line_above.setPen(gap_pen);          #line_below.setPen(gap_pen)

        ### drawing the colored rectangle. We do this for every gap except the last one
        if gap_index!= len(self.filtered_gaps)-1:    
          region_length  = self.filtered_gaps[gap_index+1][0] - gap_relative_pos #nucleotides
          region_length += self.filtered_gaps[gap_index+1][2] - self.filtered_gaps[gap_index][2] - self.filtered_gaps[gap_index][1]      #adding the small insertions in this cds; this is equal to the differences of the  total number of gaps recorded in the next gap and in this one, minus the length of this gap.          
          x_start=(total_gap_length + gap_relative_pos + gap_length  )*self.width_per_position
          width  = self.width_per_position*region_length -1    # -1 to avoid overlap 
          x_end  = x_start + width
          # filled rectangle first, then one line above and one below
          colored_rect = QtGui.QGraphicsRectItem(x_start, 0, width, self.height, self.item)   
          colored_rect.setPen(rect_pen)
          colored_rect.setBrush(rect_brush)
          line_above=QtGui.QGraphicsLineItem( x_start, 0,           x_end, 0,           self.item) #parent=
          line_below=QtGui.QGraphicsLineItem( x_start, self.height, x_end, self.height, self.item) #parent=
          line_above.setPen(outline_pen);          line_below.setPen(outline_pen)


      self.draw_text_fields()
      self.draw_annotations()

      return  

    def draw_annotations(self):
      for annotation_name in self.annotations:
#        print "ANNOTATION   ", annotation_name, self.annotations[annotation_name]
        if self.annotations[annotation_name]:     ## to be able to turn off default tracks we add this control
          if self.colors.has_key(annotation_name):  color=self.colors[annotation_name]
          else:                                     color='#000000'
          pen_this_annotation=    QtGui.QPen(     QtGui.QColor(color) )
          brush_this_annotation=  QtGui.QBrush(   QtGui.QColor(color) )
          for ann_index, ann_start in enumerate(self.annotations[annotation_name]['start']):
            ann_end= self.annotations[annotation_name]['end'][ann_index]
            x_start= (ann_start)*self.width_per_position
            x_end=   (ann_end)  *self.width_per_position

            y0_relative=self.annotations[annotation_name]['y0'][ann_index]
            y1_relative=self.annotations[annotation_name]['y1'][ann_index]            
            y0_absolute=self.height*(1.0 - y0_relative)      #turning upside down as qt4 likes
            y1_absolute=self.height*(1.0 - y1_relative)      #turning upside down as qt4 likes
            
            y_min  =y1_absolute;    y_max=y0_absolute

            if x_start != x_end:
              ## drawing rectangle
              colored_rect = QtGui.QGraphicsRectItem(x_start, y_min, x_end-x_start    , y_max-y_min, self.item)   ## no border
              colored_rect.setPen(pen_this_annotation)
              colored_rect.setBrush(brush_this_annotation)
            else: 
              ## drawing line
              colored_line= QtGui.QGraphicsLineItem( x_start, y_min, x_end, y_max, self.item)
              colored_line.setPen(pen_this_annotation)


class SmallExonViewFace(ExonViewFace):
    default_attributes={
      'order_text':['id', 'chromosome', 'boundaries', 'text'],
      'height':20, 'width':100, 'width_per_position':None, 'min_for_gaps':6,
      'font':'Courier New', 'font_size':12, 'shrink_font_to_fit':True, 'font_margin':4,
      'colors':{'bkg':'#4444AA', 'outline':'#000000', 'text':'#000000', 'gap':'#999999', 'introns':'#000000', 'frameshifts':'#EE2222' },
      'printed':{'id':0, 'chromosome':1, 'text':0, 'boundaries':0},
      'annotations':{'introns':True, 'frameshifts':True}
    }


#print globals()

"""
def remember():                                                
#gene_brick_height/6, offset_for_id, gene_brick_height*2/3) #parent item: bkg

      offset_for_additional_here=offset_for_additional*len(self.gene.additional_features)      


      if not opt['no_id']:
        ### rectangle for prediction id
        chrom_id_rect = QtGui.QGraphicsRectItem(0, gene_brick_height/6, offset_for_id, gene_brick_height*2/3) #parent item: bkg
        chrom_id_rect.setParentItem(self.item)
        ## text for prediction id
        obj=QtGui.QGraphicsSimpleTextItem() 
        obj.setText(self.gene.id)
        shrink_font_to_box(obj, offset_for_id, gene_brick_height*2/3, category='prediction id')  ## modifying font size in place
        

        obj.setPos( (offset_for_id-obj.boundingRect().width())/2 , (gene_brick_height-obj.boundingRect().height())/2 ) #centering the text
        obj.setParentItem(chrom_id_rect)

      for index, x in enumerate( self.gene.additional_features ):
        #write(self.gene.full_id+' drawing secis '+' '+str(offset_for_additional), 1)
        secis_rect = QtGui.QGraphicsRectItem(gene_brick_width+offset_for_id+offset_for_additional*index, gene_brick_height/8, offset_for_additional-1, gene_brick_height*6/8) #parent item: bkg
        secis_rect.setParentItem(self.item)
        secis_rect.setBrush(QtGui.QBrush(QtGui.QColor( x.color()    )))
        #pen=QtGui.QPen(); pen.width=3
        #pen.setColor(QtGui.QColor("#FFFFFF"))
        #secis_rect.setPen(QtGui.QPen(color="#FFFFFF", width=2) )
        obj=QtGui.QGraphicsSimpleTextItem()         
        obj.setText(x.text)
        obj.setBrush( QtGui.QBrush(QtGui.QColor( "#FFFFFF"  ) ))
        reduce_font_if_necessary(obj, offset_for_additional_here, gene_brick_height*6/8, category='secis description')
        obj.setParentItem(secis_rect)
        obj.setPos(  gene_brick_width+offset_for_id+offset_for_additional*index+ (offset_for_additional-obj.boundingRect().width())/2 , (gene_brick_height-obj.boundingRect().height())/2 ) #centering the text

      #drawing line to represent the 100 coverage for profile
      line_for_full_coverage= QtGui.QGraphicsLineItem(offset_for_id, 0, offset_for_id+gene_brick_width, 0)
      line_for_full_coverage.setParentItem(self.item)
      ## rectangle for gene brick
      gene_brick=QtGui.QGraphicsRectItem(offset_for_id+self.gene.relative_boundaries()[0]*gene_brick_width, 0, (self.gene.relative_boundaries()[1]-self.gene.relative_boundaries()[0])*gene_brick_width, gene_brick_height) 
      gene_brick.setBrush(QtGui.QBrush(QtGui.QColor(self.gene.color())))
      gene_brick.setParentItem(self.item)

      #drawing lines for introns
      if not opt['I'] and len(self.gene.exons)>1:
        cds_so_far=0
        tot_cds=self.gene.length()
        for exon_index in range(len(self.gene.exons[:-1])):
          start, end = self.gene.exons[exon_index]
          cds_so_far+= end-start+1
          aa_so_far = cds_so_far/3.0
          if self.gene.strand=='+':     intron_length= self.gene.exons[exon_index+1][0] - end -1   #######
          elif self.gene.strand=='-':   intron_length= start - self.gene.exons[exon_index+1][1] -1   #######
          x_intron  =  self.gene.relative_position_in_ali_of(aa_so_far) * gene_brick_width  +  offset_for_id
          line_intron= QtGui.QGraphicsLineItem( x_intron, 1, x_intron, gene_brick_height-2)
          
          color=opt['intron_color'] # '#FFFFFF' #white
          if intron_length<=5: color='#EE0000' #red for frameshifts
          line_intron.setPen(QtGui.QPen(QtGui.QColor(color)))
          line_intron.setZValue(1)
          line_intron.setParentItem(gene_brick)          

      if opt['f']:
        for feature_name in self.gene.graphical_features: 
          for aa_position in self.gene.graphical_features[feature_name]:
            x_feature= self.gene.relative_position_in_ali_of(aa_position) * gene_brick_width  +  offset_for_id
            line_feature = QtGui.QGraphicsLineItem( x_feature, 1, x_feature, gene_brick_height-2)
            if graphical_features_colors.has_key(feature_name):
              line_feature.setPen(QtGui.QPen(QtGui.QColor( graphical_features_colors[feature_name]  )))
            line_feature.setZValue(1.1)
            line_feature.setParentItem(gene_brick)
            
      if not opt['T']:
      ## text for chromosome
        obj=QtGui.QGraphicsSimpleTextItem() 
        obj.setPos(offset_for_id+1, 0)
        obj.setText(self.gene.chromosome)
        font=obj.font();  font.setPointSize(opt['bsize']);     obj.setFont(font) #setting default text size
        reduce_font_if_necessary(obj, gene_brick_width, category='text for chromosome')
        obj.setZValue(2)
        obj.setParentItem(gene_brick)
      ## text for positions
        obj=QtGui.QGraphicsSimpleTextItem() 
        obj.setPos(offset_for_id+1,gene_brick_height/2 +1 )
        obj.setText(join([str(i) for i in self.gene.boundaries()], self.gene.strand)) # e.g. 1+101 (positive strand) 45-80 (negative strand
        font=obj.font();  font.setPointSize(opt['bsize']);     obj.setFont(font) #setting default text size      
        reduce_font_if_necessary(obj, gene_brick_width, category='positions text')
        obj.setZValue(2.1)
        obj.setParentItem(gene_brick)

    def _width(self):        return self.item.rect().width()
    def _height(self):       return self.item.rect().height()

"""
