# $Id: GeneIDServices.pm,v 1.4 2005-05-11 15:20:16 gmaster Exp $
#
# INBPerl module for INB::GRIB::geneid::MobyParser
#
# This file is an instance of a template written 
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#
# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::GeneIDServices  - Package for parser the Moby message to call the geneid gene finding program .

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

Francisco Camara, fcamara@imim.es

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

package INB::GRIB::Services::GeneIDServices;

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
  &runGeneID
  &runGeneIDGFF 
);

our $VERSION = '1.0';


# Preloaded methods go here.

###############################################################################

=head2 _do_query_GeneID_CGI

 Title   : _do_query_GeneID_CGI
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_GeneID_CGI ($query);
         : 
         : donde:
         :   $query es un XML::DOM::Node que contiene el arbol xml que
         :   engloba:  
         :     <moby:mobyData queryID='1'>...</moby:mobyData>
         : 
         : Call the GeneID CGI script and return HTML output
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

sub _do_query_GeneID_CGI {
	# $query es un objeto DOM::Node con la informacion de una query biomoby 
	my $query = shift @_;
        my $MOBY_RESPONSE = "";             # set empty response
	# Aqui escribimos las variables que necesitamos para la funcion. 
	my $nucleotide;  
	my $sequenceIdentifier;

	# get the queryID attribute of the queryInput
        my $queryID = getInputID($query);  
	# Obtenemos todos los articulos que forman parte de la query. 
	my @input_articles = getArticles($query);

	# Tratamos a cada uno de los articulos
	foreach my $article(@input_articles) {       

		# El articulo es una tupla que contiene el nombre de este 
		# y su texto xml. 

		my ($articleName, $DOM) = @{$article}; # get the named article
		
		# Si le hemos puesto nombre a los articulos del servicio,  
		# podemos recoger a traves de estos nombres el valor. 

		if ($articleName eq "sequence") { 

		    # Test get the sequence identifier another way !
		    
		    my @articles = ($DOM);
		    my @ids = getSimpleArticleIDs (\@articles);
		    my $seqId = "Moby";
		     if (@ids > 0) {
			 $seqId = $ids[0];
		    }
		    
		    # Los contenidos los devuelve como una lista, dado que 
		    # el objeto de la ontologia podria tener una relacion
		    # "has" n-aria. Bien, en nuestro caso solo habia un peptido. 

		    # The Sequence as a string

		    ($nucleotide) = getNodeContentWithArticle($DOM, "String", "SequenceString");
		    # Lo que hacemos aqui es limpiar un sting de caracteres raros 
		    # (espacios, \n, ...) pq nadie asegura que no los hayan. 
		    $nucleotide =~ s/\W+//sg; # trim trailing whitespace

		}
	}
	# Una vez recogido todos los parametros necesarios, llamamos a 
	# la funcion que nos devuelve el report. 

	my $report = GeneID_call_CGI (nucleotide  => $nucleotide, seqIdentifier => $sequenceIdentifier);

	# Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX 
	# nos queda encapsularla en un Objeto bioMoby. Esta operacio 
	# la podriamos realizar en una funcion a parte si fuese compleja.  
my $input = <<PRT;
<moby:text-html namespace='' id=''>
<![CDATA["
$report
"]]>
</moby:text-html>
PRT
	# Bien!!! ya tenemos el objeto de salida del servicio , solo nos queda
	# volver a encapsularlo en un objeto biomoby de respuesta. Pero 
	# en este caso disponemos de una funcion que lo realiza. Si tuvieramos 
	# una respuesta compleja (de verdad, esta era simple ;) llamariamos 
	# a collection response. 
	# IMPORTANTE: el identificador de la respuesta ($queryID) debe ser 
	# el mismo que el de la query. 
	$MOBY_RESPONSE .= simpleResponse($input, "report", $queryID);
	
        return $MOBY_RESPONSE;
}


=head2 _do_query_GeneID

 Title   : _do_query_GeneID
         : 
         : private function (NOT EXPORTED)
         : 
 Usage   : my $query_response = _do_query_GeneID($query);
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

sub _do_query_GeneID {
    # $queryInput_DOM es un objeto DOM::Node con la informacion de una query biomoby 
    my $queryInput_DOM = shift @_;
    # $_format is the type of output that returns GeneID (e.g. GFF)
    my $_format        = shift @_;
    
    my $MOBY_RESPONSE = "";     # set empty response
    
    # Aqui escribimos las variables que necesitamos para la funcion. 
    my $profile;
    my $strands;
    
    # Variables that will be passed to GeneID_call
    my %sequences;
    my %parameters;
    
    my $queryID  = getInputID ($queryInput_DOM);
    my @articles = getArticles($queryInput_DOM);
    
    # Get the parameters
    
    ($profile) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "profile");
    if (not defined $profile) {
	# Default is "Human"
	$profile = "Human";
    }
    
    ($strands) = getNodeContentWithArticle($queryInput_DOM, "Parameter", "strands");
    if (not defined $strands) {
	# Default is running GeneID on both strands
	$strands = "Both";
    }
    
    # Add the parsed parameters in a hash table
    
    $parameters{profile} = $profile;
    $parameters{strands} = $strands;
    
    # Tratamos a cada uno de los articulos
    foreach my $article (@articles) {       
	
	# El articulo es una tupla que contiene el nombre de este 
	# y su texto xml. 
	
	my ($articleName, $DOM) = @{$article}; # get the named article
	
	# Si le hemos puesto nombre a los articulos del servicio,  
	# podemos recoger a traves de estos nombres el valor. 	
	
	if ($articleName eq "sequences") { 
	    
	    if (isSimpleArticle ($DOM)) {

		# print STDERR "sequences tag is a simple article...\n";
		
		my $sequenceIdentifier;
		my $nucleotide;
		
		my @articles = ($DOM);
		my @ids = getSimpleArticleIDs (\@articles);
		if (@ids > 0) {
		    $sequenceIdentifier = $ids[0];
		}
		else {
		    print STDERR "Error - no sequence identifier!!!\n";
		}
		
		# Los contenidos los devuelve como una lista, dado que 
		# el objeto de la ontologia podria tener una relacion
		# "has" n-aria. Bien, en nuestro caso solo habia un peptido. 
		
		# The Sequence as a string
		
		($nucleotide) = getNodeContentWithArticle($DOM, "String", "SequenceString");
		# Lo que hacemos aqui es limpiar un sting de caracteres raros 
		# (espacios, \n, ...) pq nadie asegura que no los hayan. 
		$nucleotide =~ s/\W+//sg; # trim trailing whitespace
		
		# Add the sequence into a hash table
		
		$sequences{$sequenceIdentifier} = $nucleotide;		
	    }
	    elsif (isCollectionArticle ($DOM)) {
		
		# print STDERR "sequences is a collection article...\n";
		# print STDERR "Collection DOM: " . $DOM->toString() . "\n";
		
		my @sequence_articles_DOM = getCollectedSimples ($DOM);
		
		foreach my $sequence_article_DOM (@sequence_articles_DOM) {
		    
		    my ($sequenceIdentifier) = getSimpleArticleIDs ( [ $sequence_article_DOM ] );

		    # print STDERR "Sequence DOM: " . $sequence_article_DOM->toString() . "\n";

		    my ($nucleotide) = getNodeContentWithArticle($sequence_article_DOM, "String", "SequenceString");
		    # Lo que hacemos aqui es limpiar un sting de caracteres raros 
		    # (espacios, \n, ...) pq nadie asegura que no los hayan.
		    $nucleotide =~ s/\W+//sg; # trim trailing whitespace
		    
		    if (length ($nucleotide) < 1) {
			print STDERR "nucleotide sequence not parsed properly...\n";
		    }

		    # print STDERR "sequenceIdentifier: $sequenceIdentifier\n";
		    # print STDERR "nucleotide length: " . length ($nucleotide) . "\n";
		    
		    # Add the sequence into a hash table
		    
		    $sequences{$sequenceIdentifier} = $nucleotide;
		}
	    }
	    else {
		print STDERR "It is not a simple or collection article...\n";
		print STDERR "DOM: " . $DOM->toString() . "\n";
	    }
	    
	} # End parsing sequences article tag
	
    } # Next article
    
    # Check that we have parsed properly the sequences
    
    if ((keys (%sequences)) == 0) {
	print STDERR "Error, can't parsed any sequences...\n";
    }
    
    # Una vez recogido todos los parametros necesarios, llamamos a 
    # la funcion que nos devuelve el report. 	
    
    my $report = GeneID_call (sequences  => \%sequences, format => $_format, parameters => \%parameters);
    
    # Ahora que tenemos la salida en el formato de la aplicacion XXXXXXX 
    # nos queda encapsularla en un Objeto bioMoby. Esta operacio 
    # la podriamos realizar en una funcion a parte si fuese compleja.  
    
    my $input = <<PRT;
<moby:$_format namespace='' id=''>
<![CDATA["
$report
"]]>
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




=head2 runGeneID 

 Title   : runGeneID 
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

sub runGeneID {

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

	# Para cada query ejecutaremos el _execute_query.
        foreach my $query(@queries){
		# En este punto es importante recordar que el objeto $query 
		# es un XML::DOM::Node, y que si queremos trabajar con 
		# el mensaje de texto debemos llamar a: $query->toString() 
		my $query_response = _do_query_GeneID($query);
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

=head2 runGeneIDGFF

 Title   : runGeneIDGFF 
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
         : queries en GFF formate. Es decir, un mensaje xml de la forma:
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

sub runGeneIDGFF {

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
	# The output format for this service is GFF
	#
        my $_format = "GFF";

	# Para cada query ejecutaremos el _execute_query.
        foreach my $queryInput (@queries){

	    # En este punto es importante recordar que el objeto $query 
	    # es un XML::DOM::Node, y que si queremos trabajar con 
	    # el mensaje de texto debemos llamar a: $query->toString() 
	    
	    # my $query_str = $queryInput->toString();
	    # print STDERR "query text: $query_str\n";

	    my $query_response = _do_query_GeneID ($queryInput, $_format);
	    
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
