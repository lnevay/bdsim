/*
   bison grammar for the gmad parser
   Ilya Agapov, 2005
*/

//TODO : implement loops!


%{

  extern int line_num;
  extern char* yyfilename;


  const int DEBUG = 0; // print debug info like parsing output etc.
  const int ECHO_GRAMMAR = 0; // print grammar rule expamsion (for debugging)
  const int VERBOSE = 0; // print warnings and errors
  const int VERBOSE_EXPAND = 0; // print the process of line expansion 
  const int INTERACTIVE =0; // print output of commands (like in interactive mode)

#include "parser.h"

  int execute = 1;
  int element_count = 1; // for samplers , ranges etc.

%}


%union{
  double dval;
  int ival;
  struct symtab *symp;
  char *str;
  struct Array *array;
}


%left '+' '-'
%left '*' '/'
%left '^' '<' '>' NE LE GE
%nonassoc UMINUS
%nonassoc UPLUS

%token <dval> NUMBER
%token <symp> VARIABLE VECVAR FUNC 
%token <str> STR
%token MARKER ELEMENT DRIFT RF DIPOLE SBEND QUADRUPOLE SEXTUPOLE OCTUPOLE MULTIPOLE 
%token SOLENOID COLLIMATOR RCOL ECOL LINE SEQUENCE SPOILER ABSORBER LASER TRANSFORM3D
%token VKICK HKICK KICK
%token PERIOD APERTURE FILENAME GAS PIPE TUNNEL MATERIAL
%token BEAM OPTION PRINT RANGE STOP USE VALUE ECHO PRINTF SAMPLE CSAMPLE
%token IF ELSE BEGN END LE GE NE FOR
%token CUT

%type <dval> aexpr
%type <dval> expr
%type <symp> assignment
%type <array> vecexpr;
%type <array> vect;
%type <str> use_parameters;
%type <ival> extension;
%type <symp> sample_options
%type <symp> csample_options
%type <symp> gas_options
%%

input : 
      | input stmt ';'
       { 
	 if(ECHO_GRAMMAR) printf("input -> input stmt ';' \n");
       }
/*       | error input';' */
/* { */
/*   yyerrok; */
/* } */
;


stmt :          if_clause stmt { if(ECHO_GRAMMAR) printf("stmt -> IF '(' aexpr ')' stmt\n" ); execute = 1;}
              | if_clause stmt ELSE stmt  
                { if(ECHO_GRAMMAR) printf("stmt -> IF '(' bool_expr ')' ELSE stmt \n" ); }
              | atomic_stmt  { if(ECHO_GRAMMAR) printf("stmt -> atomic_stmt \n"); }
              | BEGN input END { if(ECHO_GRAMMAR) printf("stmt -> '{' stmt ';' atomic_stmt '}' \n"); }
              
;

if_clause: IF '(' aexpr ')' {if( ($3 > 0) && (execute > 0) ) execute = 1; else execute = 0;}
;

atomic_stmt : 
            | expr { if(ECHO_GRAMMAR) printf("atomic_stmt -> expr\n"); }
            | command  { if(ECHO_GRAMMAR) printf("atomic_stmt -> command\n"); }
            | decl  { if(ECHO_GRAMMAR) printf("atomic_stmt -> decl\n"); }
            | error
              {
		//yyerror(" : some error message\n"); 
		if(ECHO_GRAMMAR) printf("\natomic_stmt -> error\n");
	      }
;



decl : VARIABLE ':' marker
       {
	 if(execute)  {
	   if(ECHO_GRAMMAR) printf("decl -> VARIABLE (%s) : marker\n",$1->name);
	   // check parameters and write into element table
	   write_table(params,$1->name,_MARKER);
	   params.flush();
	 }
       }
     | VARIABLE ':' drift
       {
	 if(execute) {
	   if(ECHO_GRAMMAR) printf("decl -> VARIABLE (%s) : drift\n",$1->name);
	   // check parameters and write into element table
	   write_table(params,$1->name,_DRIFT);
	   params.flush();
	 }
       } 
     | VARIABLE ':' rf
       {
	 if(execute) {
	   if(ECHO_GRAMMAR) printf("decl -> VARIABLE (%s) : rf\n",$1->name);
	   // check parameters and write into element table
	   write_table(params,$1->name,_RF);
	   params.flush();
	 }
       } 
     | VARIABLE ':' sbend
       {  
	 if(execute) {
	   if(ECHO_GRAMMAR) printf("decl -> VARIABLE (%s) : sbend\n",$1->name);
	   // check parameters and write into element table
	   write_table(params,$1->name,_SBEND);
	   params.flush();
	 }
       }
    | VARIABLE ':' vkick
       {  
	 if(execute) {
	   if(ECHO_GRAMMAR) printf("decl -> VARIABLE (%s) : vkick\n",$1->name);
	   // check parameters and write into element table
	   write_table(params,$1->name,_VKICK);
	   params.flush();
	 }
       }
    | VARIABLE ':' hkick
       {  
	 if(execute) {
	   if(ECHO_GRAMMAR) printf("decl -> VARIABLE (%s) : hkick\n",$1->name);
	   // check parameters and write into element table
	   write_table(params,$1->name,_HKICK);
	   params.flush();
	 }
       }
     | VARIABLE ':' quad
       {
	 if(execute)       
	   {
	     if(ECHO_GRAMMAR) printf("decl -> VARIABLE : quad %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_QUAD);
	     params.flush();
	   }
       }
     | VARIABLE ':' sextupole
       {
	 if(execute)
	   {
	     if(ECHO_GRAMMAR) printf("decl -> VARIABLE : sext %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_SEXTUPOLE);
	     params.flush();
	   }
       }
     | VARIABLE ':' octupole
       {
	 if(execute)
	   {
	     if(ECHO_GRAMMAR) printf("VARIABLE : octupole %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_OCTUPOLE);
	     params.flush();
	   }
       }
     | VARIABLE ':' rcol
       {
	 if(execute)
	   {
	     if(ECHO_GRAMMAR) printf("VARIABLE : rcol %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_RCOL);
	     params.flush();
	   }
       }
     | VARIABLE ':' ecol
       {
	 if(execute)
	   {
	     if(ECHO_GRAMMAR) printf("VARIABLE : ecol %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_ECOL);
	     params.flush();
	   }
       }
     | VARIABLE ':' multipole
       {
	 if(execute)
	   {	 
	     if(ECHO_GRAMMAR) printf("VARIABLE : multipole %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_MULT);
	     params.flush();	 
	   }
       }
     | VARIABLE ':' element
       {
	 if(execute)
	   {	 
	     if(ECHO_GRAMMAR) printf("VARIABLE : element %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_ELEMENT);
	     params.flush();	 
	   }
       }
     | VARIABLE ':' laser
       {
	 if(execute)
	   {	 
	     if(ECHO_GRAMMAR) printf("VARIABLE : laser %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_LASER);
	     params.flush();	 
	   }
       }
     | VARIABLE ':' transform3d
       {
	 if(execute)
	   {	 
	     if(ECHO_GRAMMAR) printf("VARIABLE : transform3d %s \n",$1->name);
	     // check parameters and write into element table
	     write_table(params,$1->name,_TRANSFORM3D);
	     params.flush();	 
	   }
       }
     | VARIABLE ':' line 
       {
	 if(execute)
	   {
	     // create entry in the main table and add pointer to the parsed sequence
	     if(ECHO_GRAMMAR) printf("VARIABLE : LINE %s\n",$1->name);
	     write_table(params,$1->name,_LINE,new list<struct Element>(tmp_list));
	     tmp_list.erase(tmp_list.begin(), tmp_list.end());
	   }
       }
     | VARIABLE ':' sequence
       {
	 if(execute)
	   {
	   }
       }
     | VARIABLE ':' extension
       {
	 if(execute)
	   {
	     if(ECHO_GRAMMAR) printf("decl -> VARIABLE : VARIABLE, %s  :  %s\n",$1->name, typestr($3));
	     if($3 != _NONE)
	       {
		 write_table(params,$1->name,$3);
	       }
	     params.flush();
	   }
       }
     | VARIABLE ':' material
       {
	 if(execute)
	   {
	     if(ECHO_GRAMMAR) printf("decl -> VARIABLE : Material, %s \n",$1->name);
	     write_table(params,$1->name,_MATERIAL);
	     params.flush();
	   }
       }
;

marker : MARKER ;

drift : DRIFT parameters
;

rf : RF parameters
;

sbend : SBEND parameters
;

vkick : VKICK parameters
;

hkick : HKICK parameters
;

quad : QUADRUPOLE parameters
;

sextupole : SEXTUPOLE parameters
;

octupole : OCTUPOLE parameters
;

multipole : MULTIPOLE parameters
{
  if(DEBUG) print(params);
}
;

ecol : ECOL parameters
;

rcol : RCOL parameters
;

laser : LASER parameters
;

transform3d : TRANSFORM3D parameters
;

element : ELEMENT parameters
;

material : MATERIAL parameters
;

extension : VARIABLE parameters
            {
	      if(execute)
		{	 
		  if(ECHO_GRAMMAR) printf("extension : VARIABLE parameters   -- %s \n",$1->name);
		  list<struct Element>::iterator it = element_lookup($1->name);
		  if(it == NULL)
		    {
		      if(VERBOSE) printf("type %s has not been defined\n",$1->name);
		      $$ = _NONE;
		    }
		  else
		    {
		      // inherit properties from the base type
		      $$ = (*it).type;
		      inherit_properties(*it);
		    }
		  
		}
	    }
;

parameters: 
          | parameters ',' VARIABLE '=' aexpr
            {
	      if(execute)
		{
		  if(DEBUG) printf("parameters, VARIABLE(%s) = aexpr(%.10g)\n",$3->name,$5);
		  if(!strcmp($3->name,"l")) { params.l = $5; params.lset = 1;} 
		    else
		  if(!strcmp($3->name,"k0")) { params.k0 = $5; params.k0set = 1;}
		    else 
		  if(!strcmp($3->name,"k1")) { params.k1 = $5; params.k1set = 1;} 
		    else
		  if(!strcmp($3->name,"k2")) { params.k2 = $5; params.k2set = 1;}
		    else 
		  if(!strcmp($3->name,"k3")) { params.k3 = $5; params.k3set = 1;}
		    else 
		  if(!strcmp($3->name,"angle")) { params.angle = $5; params.angleset = 1;}
		    else
		  if(!strcmp($3->name,"aper") ||!strcmp($3->name,"aperture") ) 
			      { params.aper = $5; params.aperset = 1;}
		    else
		  if(!strcmp($3->name,"outR") ) { params.outR = $5; params.outRset = 1;}
		  if(!strcmp($3->name,"xsize") ) { params.xsize = $5; params.xsizeset = 1;}
		    else
		  if(!strcmp($3->name,"ysize") ) { params.ysize = $5; params.ysizeset = 1;}
		    else
		  if(!strcmp($3->name,"tilt")) { params.tilt = $5; params.tiltset = 1;}
		    else
		  if(!strcmp($3->name,"x")) {params.xdir = $5; params.xdirset = 1;} // x direction
		    else
		  if(!strcmp($3->name,"y")) {params.ydir = $5; params.ydirset = 1;} // y direction 
		    else
		  if(!strcmp($3->name,"z")) {params.zdir = $5; params.zdirset = 1;} // z direction 
		    else
		  if(!strcmp($3->name,"phi")) {params.phi = $5; params.phiset = 1;}  // polar angle
		    else
		  if(!strcmp($3->name,"theta"))  {params.theta = $5; params.thetaset = 1;} 
		  // azimuthal angle
		    else
		  if(!strcmp($3->name,"psi"))  {params.psi = $5; params.psiset = 1;} // 3rd  angle
		  else
		  if(!strcmp($3->name,"ez"))  {params.ez = $5; params.ezset = 1;} // rf voltage
		    else
		  if(!strcmp($3->name,"fint")) {;} // fringe field parameters
		    else
		  if(!strcmp($3->name,"fintx")) {;}  //
                    else
		  if(!strcmp($3->name,"e1")) {;}  //
                    else
		  if(!strcmp($3->name,"e2")) {;}  //
                    else
		  if(!strcmp($3->name,"hgap")) {;}  //
		    else
		  if(!strcmp($3->name,"density")) {;}  //
		    else
		  if(!strcmp($3->name,"A")) {params.A = $5; params.Aset = 1;}  //
		    else
		  if(!strcmp($3->name,"Z")) {params.Z = $5; params.Zset = 1;}  //
		    else
		  if(!strcmp($3->name,"density")) {params.density = $5; params.densityset = 1;}  //
                    else
		  if(!strcmp($3->name,"T")) {params.temper = $5; params.temperset = 1;}  //
		    else
		  if(VERBOSE) printf("Warning : unknown parameter %s\n",$3->name);
		  
		}
	    }
           | parameters ',' VARIABLE '=' vecexpr
             {
	       if(execute) 
		 {
		   if(DEBUG) printf("params,VARIABLE (%s) = vecexpr (%d)\n",$3->name,$5->size);
		   if(!strcmp($3->name,"knl")) 
		     {
		       params.knlset = 1;
		       set_vector(params.knl,$5);
		       delete[] $5->data;
		     } 
		   else
		     if(!strcmp($3->name,"ksl")) 
		       {
			 params.kslset = 1;
			 set_vector(params.ksl,$5);
			 delete[] $5->data;
		       }
		     else 	  
		       if(VERBOSE) printf("unknown parameter %s\n",$3->name);
		 }
	     }         
           | VARIABLE '=' vecexpr
             {
	       if(execute) 
		 {
		   if(DEBUG) printf("VARIABLE (%s) = vecexpr\n",$1->name);
		   if(!strcmp($1->name,"knl")) 
		     {
		       params.knlset = 1;
		       set_vector(params.knl,$3);
		       delete[] $3->data;
		     } 
		   else
		     if(!strcmp($1->name,"ksl")) 
		       {
			 params.kslset = 1;
			 set_vector(params.ksl,$3);
			 delete[] $3->data;
		       }
		     else 	  
		       if(VERBOSE) printf("unknown parameter %s\n",$1->name);
		 }         
	     }
          | VARIABLE '=' aexpr
            {
	      if(execute)
		{
		  if(DEBUG) printf("VARIABLE (%s) = aexpr(%.10g)\n",$1->name,$3);
		}
	    }
          | parameters ',' VARIABLE '=' STR
             {
	       if(execute) 
		 {
		   if(DEBUG) printf("params,VARIABLE (%s) = str (%s)\n",$3->name,$5);
		   if(!strcmp($3->name,"geometry")) 
		     {
		       params.geomset = 1;
		       strcpy(params.geometry, $5);

		     } 
		   else
		     if(!strcmp($3->name,"bmap")) 
		       {
			 params.geomset = 1;
			 strcpy(params.bmap, $5);
		       }
		   else 
		     if(!strcmp($3->name,"type")) 
		       {
			 //ignore the "type attribute for the moment"
		       }
		   else
		   if(!strcmp($3->name,"material")) 
		       {
			 params.materialset = 1;
			 strcpy(params.material, $5);
		       }
		     else 
		     	  
		       if(VERBOSE) printf("unknown parameter %s\n",$3->name);
		 }
	     }         
           | VARIABLE '=' STR
             {
	       if(execute) 
		 {
		   if(DEBUG) printf("VARIABLE (%s) = str\n",$1->name);
		   if(!strcmp($1->name,"geometry")) 
		     {
		       params.geomset = 1;
		       strcpy(params.geometry, $3);
		       
		     } 
		   else
		     if(!strcmp($1->name,"bmap")) 
		       {
			 params.geomset = 1;
			 strcpy(params.bmap, $3);
		       }
		     else 
		     if(!strcmp($1->name,"type")) 
		       {
			 //ignore the "type attribute for the moment"
		       }
		   else
		     if(!strcmp($1->name,"material")) 
		       {	 
			 params.materialset = 1;
			 strcpy(params.material, $3);
		       }
		     else 
		       if(VERBOSE) printf("unknown parameter %s\n",$1->name);
		 }         
	     }

line : LINE '=' '(' element_seq ')'           
;

sequence : SEQUENCE;


element_seq : 
            | element_seq ',' VARIABLE 
              {
		if(execute)
		  {
		    if(DEBUG) printf("matched sequence element, %s\n",$3->name);
		    // add to temporary element sequence
		    {
		      struct Element e;
		      e.name = $3->name;
		      e.type = _LINE;
		      e.lst = NULL;
		      tmp_list.push_back(e);
		    }
		  }
	      }
            | element_seq ',' VARIABLE '*' NUMBER 
              {
		if(execute)
		  {
		    if(DEBUG) printf("matched sequence element, %s * %d \n",$3->name,(int)$5);
		    // add to temporary element sequence
		    {
		      struct Element e;
		      e.name = $3->name;
		      e.type = _LINE;
		      e.lst = NULL;
		      for(int i=0;i<(int)$5;i++)
			tmp_list.push_back(e);
		    }
		  }
	      }
            | element_seq ',' NUMBER '*' VARIABLE 
              {
		if(execute)
		  {
		    if(DEBUG) printf("matched sequence element, %s * %d \n",$5->name,(int)$3);
		    // add to temporary element sequence
		    {
		      struct Element e;
		      e.name = $5->name;
		      e.type = _LINE;
		      e.lst = NULL;
		      for(int i=0;i<(int)$3;i++)
			tmp_list.push_back(e);
		    }
		  }
	      }
            | VARIABLE
              {
		if(execute)
		  {
		    if(DEBUG) printf("matched last sequence element, %s\n",$1->name);
		    // add to temporary element sequence
		    {
		      struct Element e;
		      e.name = $1->name;
		      e.type = _LINE;
		      e.lst = NULL;
		      tmp_list.push_back(e);
		    }
		  }
	      }
           | VARIABLE '*' NUMBER
              {
		if(execute)
		  {
		    if(DEBUG) printf("matched last sequence element, %s * %d\n",$1->name,(int)$3);
		    // add to temporary element sequence
		    {
		      struct Element e;
		      e.name = $1->name;
		      e.type = _LINE;
		      e.lst = NULL;
		      for(int i=0;i<(int)$3;i++)
			tmp_list.push_back(e);
		    }
		  }
	      }
            | NUMBER '*' VARIABLE
              {
		if(execute)
		  {
		    if(DEBUG) printf("matched last sequence element, %s * %d\n",$3->name,(int)$1);
		    // add to temporary element sequence
		    {
		      struct Element e;
		      e.name = $3->name;
		      e.type = _LINE;
		      e.lst = NULL;
		      for(int i=0;i<(int)$1;i++)
			tmp_list.push_back(e);
		    }
		  }
	      }
;


expr : aexpr 
       { // check type ??
	 if(ECHO_GRAMMAR) printf("expr -> aexpr\n");
	 if(execute) 
	   {
	     if(INTERACTIVE) printf ("\t%.10g\n", $1); $$=$1;
	   }
       }
     | vecexpr 
       {
	 if(ECHO_GRAMMAR) printf("expr -> vecexpr\n");
	 if(execute)
	   {
	     if(INTERACTIVE)
	       for(int i=0;i<$1->size;i++)
		 {
		   printf(" %.10g ",$1->data[i]);
		 }
	     $$ = 0;
	   } 
       }
     | assignment 
       { // check type
	 if(ECHO_GRAMMAR) printf("expr -> assignment\n");
	 if(execute)
	   {
	     if(INTERACTIVE) {
	       if($1->type == _ARRAY)
		 {
		   for(list<double>::iterator it = $1->array.begin();
		       it!=$1->array.end();it++)
		     printf ("\t%.10g", (*it));
		   printf("\n");
		 }
	       else
		 printf ("\t%.10g\n", $1->value);
	     } 
	     $$=0;
	   }
       }
;

aexpr :  NUMBER               { $$ = $1;                         }
       | VARIABLE             
         { 
	   //check type ??
	   $$ = $1->value;        
          } 
       | FUNC '(' aexpr ')'   { $$ = (*($1->funcptr))($3);       } 
       | aexpr '+' aexpr      { $$ = $1 + $3;                    }
       | aexpr '-' aexpr      { $$ = $1 - $3;                    }  
       | aexpr '*' aexpr      { $$ = $1 * $3;                    }
       | aexpr '/' aexpr      { $$ = $1 / $3;                    }
       | aexpr '^' aexpr      { $$ = pow($1,$3);                 }
       | '-' aexpr  %prec UMINUS { $$ = -$2; }
       | '+' aexpr  %prec UPLUS { $$ = $2; }
       | '(' aexpr ')'         { $$ = $2;                         }
       | '<' vecexpr ',' vecexpr '>' // scalar product
         {
	   if($2->size == $4->size)
	     {
	       $$ = 0;
	       for(int i=0;i<$2->size;i++)
		 $$ += $2->data[i] * $4->data[i];
	     }
	   else
	     {
	       if(VERBOSE) printf("vector dimensions do not match");
	       $$ = _undefined;
	     }
         } 
       // boolean stuff
        | aexpr '<' aexpr { $$ = ($1 < $3 )? 1 : 0; } 
        | aexpr LE aexpr { $$ = ($1 <= $3 )? 1 : 0; } 
        | aexpr '>' aexpr { $$ = ($1 > $3 )? 1 : 0; } 
        | aexpr GE aexpr { $$ = ($1 >= $3 )? 1 : 0; } 
        | aexpr NE aexpr { $$ = ($1 != $3 )? 1 : 0; } 
        | VARIABLE '[' VARIABLE ']' 
          { 
	    if(ECHO_GRAMMAR) printf("aexpr-> %s [ %s ]\n ",$1->name, $3->name); 
	    $$ = property_lookup($1->name,$3->name);
	  }// element attributes
 ; 

assignment :  VARIABLE '=' aexpr  
              {
		if(ECHO_GRAMMAR) printf("%s \n",$1->name);
		if(execute)
		  {
		    if($1->is_reserved)
		      printf("%s is reserved",$1->name);
		    else
		      {
			$1->value = $3; $$=$1;       
		      }
		  }
	      }
           |  VARIABLE '=' vecexpr
              {
		if(execute)
		  {
		    $1->array.erase($1->array.begin(),$1->array.end());
		    for(int i=0;i<$3->size;i++)
		      $1->array.push_back($3->data[i]);
		    $1->type = _ARRAY;
		    $$ = $1;
		    delete[] $3->data;
		    $3->size = 0;
		  }
              }
           |  VECVAR '=' vecexpr
              {
		if(execute)
		  {
		    $1->array.erase($1->array.begin(),$1->array.end());
		    for(int i=0;i<$3->size;i++)
		      $1->array.push_back($3->data[i]);
		    $$ = $1;
		    delete[] $3->data;
		    $3->size = 0;
		  }
              }
;

vecexpr :  VECVAR  
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->data = new double[$1->array.size()];
	      $$->size = $1->array.size();
	      //array_list.push_back($$);
	      list<double>::iterator it = 0;
	      int i = 0;
	      for(it=$1->array.begin();it!=$1->array.end();it++)
		{
		  $$->data[i++] = (*it);
		}
	    }
        } 
        | vect
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->data = new double[$1->size];
	      $$->size = $1->size;
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$1->size;i++)
		{
		  $$->data[i] = $1->data[i];
		}
	      
	      
	      // erase data in vect
	      
	      delete[] $1->data;
	      $1->size = 0;
	    }
	}
       | vecexpr '+' vecexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = ($1->size < $3->size )? $1->size : $3->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $1->data[i] + $3->data[i];
		}
	      
	      
	      // erase data in vect
	      
	      delete[] $1->data;
	      delete[] $3->data;
	      $1->size = 0;
	      $3->size = 0;
	    }
        }
      | vecexpr '-' vecexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = ($1->size < $3->size )? $1->size : $3->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $1->data[i] - $3->data[i];
		}
	      
	      
	      // erase data in vect
	      
	      delete[] $1->data;
	      delete[] $3->data;
	      $1->size = 0;
	      $3->size = 0;
	    }
	}
       | vecexpr '+' aexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = $1->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $1->data[i] + $3;
		}
	      
	      // erase data in vect
	      
	      delete[] $1->data;
	      $1->size = 0;
	    }
	}

      | vecexpr '*' aexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = $1->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $1->data[i] * $3;
		}
	      
	      // erase data in vect
	      
	      delete[] $1->data;
	      $1->size = 0;
	    }
	}
      | vecexpr '/' aexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = $1->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $1->data[i] / $3;
		}
	      
	      // erase data in vect
	      
	      delete[] $1->data;
	      $1->size = 0;
	    }
	}
       | aexpr '+' vecexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = $3->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $3->data[i] + $1;
		}
	      
	      // erase data in vect
	      
	      delete[] $3->data;
	      $3->size = 0;
	    }
	}
       | aexpr '-' vecexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = $3->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $1 - $3->data[i];
		}
	      
	      // erase data in vect
	      
	      delete[] $3->data;
	      $3->size = 0;
	    }
	}
      | aexpr '*' vecexpr
        {
	  if(execute)
	    {
	      $$ = new struct Array;
	      $$->size = $3->size;
	      $$->data = new double[$$->size];
	      //array_list.push_back($$);
	      
	      for(int i=0;i<$$->size;i++)
		{
		  $$->data[i] = $1 * $3->data[i];
		}
	      
	      // erase data in vect
	      
	      delete[] $3->data;
	      $3->size = 0;
	    }
	}

;

vect : '{' numbers '}' 
{
  if(execute)
    {
      //printf("matched vector of size %d\n",_tmparray.size());
      $$ = new struct Array;
      $$->data = new double[_tmparray.size()];
      $$->size = _tmparray.size();
      
      //array_list.push_back(a);
      
      list<double>::iterator it;
      
      int i=0;
      for(it=_tmparray.begin();it!=_tmparray.end();it++)
	{
	  $$->data[i++] = (*it);
	}
      
      _tmparray.erase(_tmparray.begin(),_tmparray.end());
    }
}
;

numbers : 
        | numbers ',' aexpr
          {
	    if(execute)
	      _tmparray.push_back($3);
          } 
       | aexpr
         {
	   if(execute)
	     _tmparray.push_back($1);
        }
;

command : STOP             { if(execute) quit(); }
        | BEAM ',' beam_parameters
        | PRINT            { if(execute) print( element_list ); }
        | PRINT ',' LINE   { if(execute) print( beamline_list); }
        | PRINT ',' OPTION { if(execute) print(options); }
        | PRINT ',' VARIABLE 
          {
	    if(execute)
	      {
		printf("\t");
		printf("\t%.10g\n",$3->value);
	      }
	  } 
        | PRINT ',' VECVAR 
          {
	    if(execute)
	      {
		printf("\t");
		
		list<double>::iterator it;
		for(it=$3->array.begin();it!=$3->array.end();it++)
		  {
		    printf("  %.10g ",(*it));
		  }
		
		printf("\n");
	      } 
	  }
        | USE ',' use_parameters { if(execute) expand_line(current_line,current_start, current_end);}
        | OPTION  ',' option_parameters
	| ECHO STR { if(execute) printf("%s\n",$2); }
        | SAMPLE ',' sample_options 
          {
	    if(execute)
	      {  
		if(ECHO_GRAMMAR) printf("command -> SAMPLE\n");
		add_sampler("sampler",$3->name, element_count);
		element_count = 1;
		params.flush();
	      }
          }
        | CSAMPLE ',' csample_options // cylindrical sampler
          {
	    if(execute)
	      {  
		if(ECHO_GRAMMAR) printf("command -> CSAMPLE\n");
		add_csampler("sampler",$3->name, element_count,params.l, params.r);
		element_count = 1;
		params.flush();
	      }
          }
        | GAS ',' gas_options // beampipe gas
          {
	    if(execute)
	      {  
		if(ECHO_GRAMMAR) printf("command -> GAS\n");
		add_gas("gas",$3->name, element_count, params.material);
		element_count = 1;
		params.flush();
	      }
          }
        | TUNNEL ',' parameters // beampipe gas
          {
	    if(execute)
	      {  
		if(ECHO_GRAMMAR) printf("command -> TUNNEL\n");
		write_table(params,"tunnel",_TUNNEL);
		params.flush();
	      }
          }

//| PRINTF '(' fmt ')' { if(execute) printf($3,$5); }
;

use_parameters :  PERIOD '=' VARIABLE
                  {
		    if(execute)
		      {
			$$ = $3->name;
			current_line = $3->name;
			current_start = NULL;
			current_end = NULL;
		      }
                  }
                | PERIOD '=' VARIABLE ',' RANGE '=' VARIABLE '/' VARIABLE
                  {
		    if(execute)
		      {
			$$ = $3->name;
			current_line = $3->name;
			current_start = $7->name;
			current_end = $9->name;
		      }
		  }
;



sample_options: RANGE '=' VARIABLE
                { if(ECHO_GRAMMAR)  printf("sample_opt : RANGE =  %s \n",$3->name);
		  {
		    if(execute) $$ = $3;
		  }
                }
              | RANGE '=' VARIABLE '[' NUMBER ']'
                {
                  if(ECHO_GRAMMAR) printf("sample_opt : RANGE =  %s [%.10g] \n",$3->name,$5);
		    {
		      if(execute) { $$ = $3; element_count = (int)$5; }
		    }
                }
;



csample_options : VARIABLE '=' aexpr
                  {
		    if(ECHO_GRAMMAR) printf("csample_opt ->csopt , %s =  %.10g \n",$1->name,$3);
		    
		    if(execute)
		      {
			if( !strcmp($1->name,"r") ) params.r = $3;
			else if (!strcmp($1->name,"l") ) params.l = $3;
			else if(VERBOSE) 
			  printf("Warning : CSAMPLER: unknown parameter %s \n",$1->name);
		      }
		  }   
                | VARIABLE '=' STR
                  {
		    if(ECHO_GRAMMAR) printf("csample_opt -> %s =  %s \n",$1->name,$3);
		    if(execute)
		      {
			;//set_value($1->name,string($3));
		      }
		  }   
                | csample_options ',' VARIABLE '=' aexpr
                  {
		    if(ECHO_GRAMMAR) printf("csample_opt ->csopt , %s =  %.10g \n",$3->name,$5);
		    
		    if(execute)
		      {
			if( !strcmp($3->name,"r") ) params.r = $5;
			else if (!strcmp($3->name,"l") ) params.l = $5;
			else if (!strcmp($3->name,"l") ) params.l = $5;
			else if(VERBOSE) 
			  printf("Warning : CSAMPLER: unknown parameter %s at line\n",$3->name);
		      }

		  }   
                | csample_options ',' VARIABLE '=' STR
                  {
		    if(ECHO_GRAMMAR) printf("csample_opt -> %s =  %s \n",$3->name,$5);
		    if(execute)
		      ;//set_value($1->name,string($3));
		  }   
                | sample_options ',' csample_options
                  {
		    if(ECHO_GRAMMAR) printf("csample_opt -> sopt, csopt\n");
		    $$ = $1;
		  }
                | sample_options
                  {
		    if(ECHO_GRAMMAR) printf("csample_opt -> sopt\n");
		    $$ = $1;
                  }
;

gas_options : VARIABLE '=' aexpr
                  {
		    if(ECHO_GRAMMAR) printf("gas_opt -> , %s =  %.10g \n",$1->name,$3);
		    
		    if(execute)
		      {
			if( !strcmp($1->name,"r") ) params.r = $3;
			else if (!strcmp($1->name,"l") ) params.l = $3;
			else if(VERBOSE) 
			  printf("Warning : GAS: unknown parameter %s \n",$1->name);
		      }
		  }   
                | VARIABLE '=' STR
                  {
		    if(ECHO_GRAMMAR) printf("gas_opt -> %s =  %s \n",$1->name,$3);
		    if(execute)
		      {
			if( !strcmp($1->name,"material") ) 
			  {
			    strcpy(params.material ,$3);
			    params.materialset = 1;
			  }
			//set_value($1->name,string($3));
		      }
		  }   
                | gas_options ',' VARIABLE '=' aexpr
                  {
		    if(ECHO_GRAMMAR) printf("gas_opt -> , %s =  %.10g \n",$3->name,$5);
		    
		    if(execute)
		      {
			if( !strcmp($3->name,"r") ) params.r = $5;
			else if (!strcmp($3->name,"l") ) params.l = $5;
			else if(VERBOSE) 
			  printf("Warning : GAS: unknown parameter %s at line\n",$3->name);
		      }

		  }   
                | gas_options ',' VARIABLE '=' STR
                  {
		    if(ECHO_GRAMMAR) printf("csample_opt -> %s =  %s \n",$3->name,$5);
		    if(execute)
		      {
			  if( !strcmp($3->name,"material") ) 
			    {
			      strcpy(params.material ,$5);
			      params.materialset = 1;
			    }
		      }
		  }   
                | gas_options ',' RANGE '='  VARIABLE '/' VARIABLE
                  {
		    if(ECHO_GRAMMAR) printf("gas_opt -> range, csopt\n");

		  }
                | RANGE '='  VARIABLE '/' VARIABLE
                  {
		    if(ECHO_GRAMMAR) printf("gas_opt -> range\n");

                  }
                | gas_options ',' RANGE '='  VARIABLE
                  {
		    if(ECHO_GRAMMAR) printf("gas_opt -> range\n");
		    $$ = $5;
		  }
                | RANGE '='  VARIABLE
                  {
		    if(ECHO_GRAMMAR) printf("gas_opt -> range\n");
		    $$ = $3;
                  }
;

option_parameters : 
                  | option_parameters ',' VARIABLE '=' aexpr
                    {
		      if(execute)
			set_value($3->name,$5);
		    }   
                  | VARIABLE '=' aexpr
                    {
		      if(execute)
			set_value($1->name,$3);
		    } 
                  | option_parameters ',' VARIABLE '=' STR
                    {
		      if(execute)
			set_value($3->name,string($5));
		    }   
                  | VARIABLE '=' STR
                    {
		      if(execute)
			set_value($1->name,string($3));
		    }   
;

beam_parameters :
                | beam_parameters ',' VARIABLE '=' aexpr
                  {
		    if(execute)
		      set_value($3->name,$5);
		  }   
                | VARIABLE '=' aexpr
                  {
		    if(execute)
		      set_value($1->name,$3);
		  }   
                | beam_parameters ',' VARIABLE '=' STR
                  {
		    if(execute)
		      set_value($3->name,string($5));
		  }   
                | VARIABLE '=' STR
                  {
		    if(execute)
		      set_value($1->name,string($3));
		  }   
;

%%



int yyerror(char *s)
{
  printf("%s at line %d , file %s\n",s, line_num, yyfilename);
  exit(1);
}

int yywrap()
{
	return 1;
}



