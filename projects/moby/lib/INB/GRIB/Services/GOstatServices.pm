# $Id: GOstatServices.pm,v 1.1 2005-05-26 16:28:03 gmaster Exp $
#
# INBPerl module for INB::GRIB::geneid::MobyParser
#
# This file is an instance of a template written 
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::GOStatServices  - Package for parser the Moby message to call the GOstat program .

=head1 SYNOPSIS

Con este package podremos parsear las llamadas xml de BioMoby para
poder llamar al servicio runGeneid. Una vez que tengamos la salida de llamando
a las funciones de Factory.pm, podremos encapsularla a objectos BioMoby.

  # 
  # En este m�dulo parsearemos los mensajes de BioMoby para realizar la
  # llamada al programa geneid. 
  # Este modulo requiere las librerias de BioMoby. Para ver su 
  # funcionamiento en modo linea de comandos (antes de activar el servicio)
  # enviaremos directamente un mensaje <xml> en BioMoby. Por ejemplo:
  $in = <<EOF
 <?xml version='1.0' encoding='UTF-8'?>
 <moby:MOBY xmlns:moby='http://www.biomoby.org/moby-s'>
 <moby:mobyContent>
   <moby:mobyData queryID='1'>
     <moby:Simple moby:articleName='seq'>
     <moby:DNASequence namespace='Global_Keyword' id=''>
       <moby:Integer namespace="" id="" articleName="Length">126</moby:Integer>
       <moby:String namespace="" id=""  articleName="SequenceString">
    ACTGCATGCTAAAGGTACATGACCGATCGGACTGTGACTGAACCTGCATTGA 
    </moby:String>
     </moby:DNASequence>
     </moby:Simple>
     <moby:Parameter moby:articleName='arg2'><Value>swissprot</Value>
     </moby:Parameter>
   </moby:mobyData>
 </moby:mobyContent>
 </moby:MOBY>
 EOF 
  my $result = GeneID("call", $in); 

 
  # Esta llamada/s nos devuelve una variable que contiene el texto con la
  # salida del programa Geneid encapsulado en objetos Moby. 

=head1 DESCRIPTION

Este package sirve para parsear las llamadas BioMoby de entrada y salida.  
De esta forma hace de puente entre el cliente que llama el servicio y la 
aplicacion geneid.

Tipicamente la libreria que necesitamos en este m�dulo es la MOBY::CommonSubs,
que podremos ver tecleando en el terminal:
> perldoc MOBY::CommonSubs

No obstante, en esta libreria puede que no encontremos todo lo que necesitamos.
Si es as� tendiamos que utilizar un parser de XML y una libreria que lo
interprete. Estas librerias podrian ser:

use XML::DOM;
use XML::Parser

y hacer cosas como:

my $parser = new XML::DOM::Parser;
my $doc = $parser->parse($message);
my $nodes = $doc->getElementsByTagName ("moby:Simple");
my $xml   = $nodes->item(0)->getFirstChild->toString;

Atencion: recomiendo mirarse esta libreria, solo cuando CommonSubs no nos de
nuestra solucion. La raz�n de ello es que, si cambia el estandard bioMoby, en
principio las llamadas a CommonSubs no tendr�an que cambiar, pero si no las 
hacemos servir, puede que debamos modificar el c�digo.


=head1 AUTHOR

Arnaud Kerhornou, akerhornou@imim.es

=head1 COPYRIGHT

Copyright (c) 2004, Roman Roset Mayals and INB - Nodo Computacional UPC/CIRI.
 All Rights Reserved.

This module is free software; you can redistribute it and/or modify
it under the same terms as Perl itself.


=head1 DISCLAIMER

This software is provided "as is" without warranty of any kind.

=head1 APPENDIX

The rest of the documentation details each of the object
methods. Internal methods are usually preceded with a _

=cut

# Let the code begin...

package INB::GRIB::Services::GOstatServices;

use strict;
use warnings;
use Carp;

use INB::GRIB::Services::Factory; 
use MOBY::CommonSubs qw(:all);

use Data::Dumper;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# This allows declaration	use INB::UPC::Blast ':all';
# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

# Aqui pondremos las funciones que reciben el mensaje BioMoby. Ejemplo 
# 
# @EXPORT = qw( &func1 &func2);
# 
our @EXPORT = qw(
  &runGOstat
);

our $VERSION = '1.0';


# Preloaded methods go here.

###############################################################################

=head2 _do_query_GOstat

 Title   : _do_query_GOstat
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_GOstat($query);
         : 
         : donde:
         :   $query es un XML::DOM::Node que contiene el arbol xml que
         :   engloba:  
         :     <moby:mobyData queryID='1'>...</moby:mobyData>
         : 
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n
         : para una sola query.
         : Un ejemplo ser�a: 
         : 
         : <moby:mobyData moby:queryID='1'>
         :   <moby:Simple moby:articleName='report'>
         :     <moby:text-plain namespace='Global_Keyword' id=''>
	 :    ....
         :     </moby:text-plain>
         :   </moby:Simple>
         : </moby:mobyData>

=cut

sub _do_query_GOstat {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
    my $queryInput_DOM = shift @_;
    # $_format is the type of output that returns GeneID (e.g. GFF)
    my $_format        = shift @_;
    
    my $MOBY_RESPONSE = "";     # set empty response
    
    # Aqui escribimos las variables que necesitamos para la funcion. 
    # ...

    # Variables that will be passed to GeneID_call
    my @regulated_genes;
    my @array_genes;
    my %parameters;

    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);

    # Get the parameters
    
    # ...
    
    # Add the parsed parameters in a hash table
    
    # ...
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {       
	
	# El articulo es una tupla que contiene el nombre de este 
	# y su texto xml. 
	
	my ($articleName, $DOM) = @{$article}; # get the named article
	
	# Si le hemos puesto nombre a los articulos del servicio,  
	# podemos recoger a traves de estos nombres el valor.
	# Sino sabemos que es el input articulo porque es un simple/collection articulo

	# It's not very nice but taverna doesn't set up easily article name for input data so we let the users not setting up the article name of the input (which should be 'sequences')
	# In case of GeneID, it doesn't really matter as there is only one input anyway
	
	if ($articleName eq "regulated genes") { 
	    
	    if (isSimpleArticle ($DOM)) {

		print STDERR "\"regulated genes\" tag is a simple article...\n";
		
		my @articles = ($DOM);
		@regulated_genes = getSimpleArticleIDs (\@articles);

		if ((@regulated_genes == 0) || (((@regulated_genes == 1) && (not defined $regulated_genes[0])))) {
		    print STDERR "Error, no gene identifier found for the regulated genes list!\n";
		}

	    }
	    elsif (isCollectionArticle ($DOM)) {
		
		print STDERR "\"array genes\" is a collection article...\n";
		# print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		
		my @genes_articles_DOM = getCollectedSimples ($DOM);
		@regulated_genes = getSimpleArticleIDs (\@genes_articles_DOM);

		print STDERR "Got a list of " . @regulated_genes . " regulated genes\n";
	    }
	    else {
		print STDERR "It is not a simple or collection article...\n";
		print STDERR "DOM: " . $DOM->toString() . "\n";
	    }
	} # End parsing regulated genes tag
	    
	if ($articleName eq "array genes") { 
	    
	    if (isSimpleArticle ($DOM)) {
		
		print STDERR "\"array genes\" tag is a simple article...\n";
		
		my @articles = ($DOM);
		@array_genes = getSimpleArticleIDs (\@articles);
		
		if ((@array_genes == 0) || (((@array_genes == 1) && (not defined $array_genes[0])))) {
		    print STDERR "Error, no gene identifier found for the array genes list!\n";
		}
		
	    }
	    elsif (isCollectionArticle ($DOM)) {
		
		print STDERR "\"array genes\" is a collection article...\n";
		# print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		
		my @genes_articles_DOM = getCollectedSimples ($DOM);
		@array_genes = getSimpleArticleIDs (\@genes_articles_DOM);
		
		print STDERR "Got a list of " . @array_genes . " array genes\n";
	    }
	    else {
		print STDERR "It is not a simple or collection article...\n";
		print STDERR "DOM: " . $DOM->toString() . "\n";
	    }
		
	} # End parsing array genes article tag
	
    } # Next article
    
    # Una vez recogido todos los parametros necesarios, llamamos a 
    # la funcion que nos devuelve el report. 	
    
    my $report = GOstat_call (regulated_genes => \@regulated_genes, array_genes => \@array_genes, format => $_format, parameters => \%parameters);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX 
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio 
    # la podriamos realizar en una funcion a parte si fuese compleja.  
    
    my $input = <<PRT;
<moby:$_format namespace='' id=''>
<String namespace='' id=''><![CDATA[
$report
]]>
</String>
</moby:$_format>
PRT
    # Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
    # volver a encapsularlo en un objeto biomoby de respuesta. Pero 
    # en este caso disponemos de una funcion que lo realiza. Si tuvieramos 
    # una respuesta compleja (de verdad, esta era simple ;) llamariamos 
    # a collection response. 
    # IMPORTANTE: el identificador de la respuesta ($queryID) debe ser 
    # el mismo que el de la query. 

    $MOBY_RESPONSE .= simpleResponse($input, "seqFeatures", $queryID);
	
    return $MOBY_RESPONSE;
}


=head2 runGOstat

 Title   : runGOstat
 Usage   : Esta funci�n est� pensada para llamarla desde un cliente SOAP. 
         : No obstante, se recomienda probarla en la misma m�quina, antes 
         : de instalar el servicio. Para ello, podemos llamarla de la 
         : siguiente forma:
         : 
         : my $result = GeneID("call", $in);
         : 
         : donde $in es texto que con el mensaje biomoby que contendr�a
         : la parte del <tag> "BODY" del mensaje soap. Es decir, un string
         : de la forma: 
         :  
         :  <?xml version='1.0' encoding='UTF-8'?>
         :   <moby:MOBY xmlns:moby='http://www.biomoby.org/moby-s'>
         :    <moby:mobyContent>
         :      <moby:mobyData queryID='1'> 
         :      ...
         :      </moby:mobyData>
         :    </moby:mobyContent>
         :   </moby:mobyContent>
         :  </moby:MOBY>
 Returns : Devuelve un string que contiene el resultado de todas las 
         : queries. Es decir, un mensaje xml de la forma:
         :
         : <?xml version='1.0' encoding='UTF-8'?>
         : <moby:MOBY xmlns:moby='http://www.biomoby.org/moby' 
         : xmlns='http://www.biomoby.org/moby'>
         :   <moby:mobyContent moby:authority='inb.lsi.upc.es'>
         :     <moby:mobyData moby:queryID='1'>
         :       ....
         :     </moby:mobyData>
         :     <moby:mobyData moby:queryID='2'>
         :       ....
         :     </moby:mobyData>
         :   </moby:mobyContent>
         :</moby:MOBY>

=cut

sub runGOstat {

	# El parametro $message es un texto xml con la peticion.
	my ($caller, $message) = @_;        # get the incoming MOBY query XML

	# Hasta el momento, no existen objetos Perl de BioMoby paralelos 
	# a la ontologia, y debemos contentarnos con trabajar directamente 
	# con objetos DOM. Por consiguiente lo primero es recolectar la 
	# lista de peticiones (queries) que tiene la peticion. 
	# 
	# En una misma llamada podemos tener mas de una peticion, y de  
	# cada servicio depende la forma de trabajar con ellas. En este 
	# caso las trataremos una a una, pero podriamos hacer Threads para 
	# tratarlas en paralelo, podemos ver si se pueden aprovechar resultados 
	# etc.. 
	my @queries = getInputs($message);  # returns XML::DOM nodes
	# 
	# Inicializamos la Respuesta a string vacio. Recordar que la respuesta
	# es una coleccion de respuestas a cada una de las consultas.
        my $MOBY_RESPONSE = "";             # set empty response

	#
	# The moby output format for this service is text-html
	# (The GeneID output format for this service is by default GFF - right now it is hardcoded)
	#

	my $_moby_output_format   = "text-formated";

	# Para cada query ejecutaremos el _execute_query.
        foreach my $query(@queries){
	    
	    # En este punto es importante recordar que el objeto $query 
	    # es un XML::DOM::Node, y que si queremos trabajar con 
	    # el mensaje de texto debemos llamar a: $query->toString() 
	    my $query_response = _do_query_GOstat ($query, $_moby_output_format);
	    
	    # $query_response es un string que contiene el codigo xml de
	    # la respuesta.  Puesto que es un codigo bien formado, podemos 
	    # encadenar sin problemas una respuesta con otra. 
	    $MOBY_RESPONSE .= $query_response;
	}
	# Una vez tenemos la coleccion de respuestas, debemos encapsularlas 
	# todas ellas con una cabecera y un final. Esto lo podemos hacer 
	# con las llamadas de la libreria Common de BioMoby. 
	return responseHeader("genome.imim.es") 
	. $MOBY_RESPONSE . responseFooter;
}


1;

__END__
