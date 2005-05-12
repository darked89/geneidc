# $Id: Factory.pm,v 1.7 2005-05-12 09:51:16 arnau Exp $
#
# INBPerl module for INB::GRIB::geneid::Factory
#
# This file is an instance of a template written 
# by Roman Roset, INB (Instituto Nacional de Bioinformatica), Spain.
#

# POD documentation - main docs before the code


=head1 NAME

INB::GRIB::Services::Factory - Package for calling geneid service.

=head1 SYNOPSIS

  # En este package tenemos las funciones necesarias para hacer las
  # llamadas al programa Geneid  La forma de utilizarlo en esta
  # versi�n es la siguiente:

 @params = ('nucleotide' => "ATCAGAGAGCGAGCAGATGCAGAGTAGAGCGAGCGAGCGCT");
 $report = INB::GRIB::Services::Factory::factory_call(@params);

  # Esta llamada nos devuelve una variable que contiene el texto con la
  # salida del programa geneid. 

=head1 DESCRIPTION

Este package sirve para hacer la llamadas al programa geneid.
Se han utilizado las libreriras de geneid, por tanto ser�n
necesarias instalarlas previamente.

Es requisito necesario para ejecutar este m�dulo, que las variables de
entorno DIR y DATADIR est�n descritas en el sistema. Dentro de
la variable $DIR se encuentra el path absoluto a la aplicacion geneid,
mientras que en $DATADIR pondremos la localizaci�n al directorio que
contiene las bases de datos.

En siguientes versiones, en este modulo se podran a�adir las diferentes
funciones para llamar a los programas que se requieran.

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

package INB::GRIB::Services::Factory;

use strict;
use warnings;

use Carp;
use CGI;
use LWP::UserAgent;
use HTTP::Request;
use HTTP::Request::Common;

use File::Temp qw/tempfile/;

# Bioperl

use Bio::SeqIO;
use Bio::Seq;

require Exporter;

our @ISA = qw(Exporter);

# Items to export into callers namespace by default. Note: do not export
# names by default without a very good reason. Use EXPORT_OK instead.
# Do not simply export all your public functions/methods/constants.

# If you do not need this, moving things directly into @EXPORT or @EXPORT_OK
# will save memory.
our %EXPORT_TAGS = ( 'all' => [ qw() ] );

our @EXPORT_OK = ( @{ $EXPORT_TAGS{'all'} } );

# Aqui pondremos las funciones que exportaremos a otros m�dulos. Ejemplo 
# 
#  @EXPORT = qw( &func1 &func2);
# 
our @EXPORT = qw(
  &GeneID_call
  &SGP2_call
);

our $VERSION = '1.00';


###############################################################################
sub _exists_env {
	my $var_name = shift;
        croak "Environment variable \$$var_name not defined.\n" 
	if (not exists $ENV{$var_name}); 
}

BEGIN {
	## in este array escribiremos las variables de entorno que le hacen
	## falta a la aplicacion para que funcione correctamente. 
	## Ejemplo:  my @needed_vars = qw(DIR DATADIR);
	my @needed_vars = qw();
	_exists_env($_) foreach @needed_vars;
}

###############################################################################

=head2 GeneID_call_CGI  

 Title   : GeneID_call_CGI 
 Usage   : $report = GeneID_call_CGI (@params);
         : 
         : ## where @params are,
         : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
         :            'arg2'  => "query 1", ## optional
         :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

=cut

sub GeneID_call_CGI {
        my %args = @_; 
        # relleno los parametros por defecto GeneID_call_CGI (nucleotide  => $nucleotide, seqIdentifier => $sequenceIdentifier);
        my $nucleotide = $args{nucleotide} || "";
	my $sequenceIdentifier = $args{seqIdentifier} || "unknown";

        # Llama a GeneID usando el interface web
        my $diagplot_url = "http://genome.imim.es/cgi-bin/GeneID_cgi/geneid_2002/geneid_2002.cgi";
	my $agent_diag   = LWP::UserAgent->new(timeout => 30000);
        my $request_diag = POST($diagplot_url,
	       "Content-type" => "form-data",
	       "Content" => [
	       "seq"     => ">$sequenceIdentifier\n".$nucleotide,
	       "profile" => "Human",
	       "engines" => "normal",
	       "strands" => "both",
	       "format"  => "gff"
        ]);
        my $result_diag = $agent_diag->request($request_diag);
	   
	return $result_diag->content;
}


=head2 GeneID_call

 Title   : GeneID_call
 Usage   : $report = GeneID_call (@params);
         : 
         : ## where @params are,
         : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
         :            'arg2'  => "query 1", ## optional
         :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

=cut

sub GeneID_call {
        my %args = @_; 
        # relleno los parametros por defecto GeneID_call (nucleotide  => $nucleotide, seqIdentifier => $sequenceIdentifier);
        my $sequences          = $args{sequences} || undef;
	my $format             = $args{format} || "";
	my $parameters         = $args{parameters} || undef;

	# Get the parameters

	my $profile = $parameters->{profile};
	my $strands = $parameters->{strands};
	
        # Llama a GeneID en local
        my $_geneid_dir  = "/home/ug/gmaster/GeneID/geneid_2002";
        my $_geneid_bin  = "bin/geneid";
        my $_geneid_args = "";
	
        if ($format eq "GFF") {
	    $_geneid_args .= "-G";
        }
	
	if ($strands eq "Both") {
	    # Default anyway
	}
	elsif ($strands eq "Forward") {
	    $_geneid_args .= "W";
	}
	elsif ($strands eq "Reverse") {
	    $_geneid_args .= "C";
	}
	
        SWITCH: {
	    if ($profile eq "Human")         { $_geneid_args .= "P $_geneid_dir/human.param"; last SWITCH; }
	    if ($profile eq "Tetraodon")     { $_geneid_args .= "P $_geneid_dir/tetraodon.param"; last SWITCH; }
	    if ($profile eq "Drosophila")    { $_geneid_args .= "P $_geneid_dir/dros.param"; last SWITCH; }
	    if ($profile eq "Celegans")      { $_geneid_args .= "P $_geneid_dir/celegans.param"; last SWITCH; }
	    if ($profile eq "Wheat")         { $_geneid_args .= "P $_geneid_dir/wheat.param"; last SWITCH; }
	    if ($profile eq "Arabidopsis")   { $_geneid_args .= "P $_geneid_dir/arabidopsis.param"; last SWITCH; }
	    if ($profile eq "Rice")          { $_geneid_args .= "P $_geneid_dir/rice.param"; last SWITCH; }
	    if ($profile eq "Plasmodium")    { $_geneid_args .= "P $_geneid_dir/plasmodium.param"; last SWITCH; }
	    if ($profile eq "Dictyostelium") { $_geneid_args .= "P $_geneid_dir/dictyostelium.param"; last SWITCH; }
	    if ($profile eq "Aspergillus")   { $_geneid_args .= "P $_geneid_dir/aspergillus.param"; last SWITCH; }
	    if ($profile eq "Neurospora")    { $_geneid_args .= "P $_geneid_dir/neurospora.param"; last SWITCH; }
	    if ($profile eq "Cryptococcus")  { $_geneid_args .= "P $_geneid_dir/cneomorfans.param"; last SWITCH; }
	    if ($profile eq "Coprinus")      { $_geneid_args .= "P $_geneid_dir/cinereus.param"; last SWITCH; }
	    # Default is Human
	    $_geneid_args .= "P $_geneid_dir/human.param";
	}
        
	# Generate a temporary file locally with the sequence(s) in FASTA format
	# locally, ie not on a NFS mounted directory, for speed sake

	my ($seq_fh, $seqfile) = tempfile("/tmp/moby/GENEID.XXXXXX", UNLINK => 0);

	# Bioperl sequence factory
	
	my $sout = Bio::SeqIO->new (
				    -fh     => $seq_fh,
				    -format => 'fasta'
				    );
	
	my @seqIds = keys (%$sequences);

	foreach my $sequenceIdentifier (@seqIds) {

	    my $nucleotide = $sequences->{$sequenceIdentifier};

	    # bioperl object
	    
	    my $seqobj = Bio::Seq->new (
					-display_id => $sequenceIdentifier,
					-seq        => $nucleotide
					);
	    
	    $sout->write_seq ($seqobj);
	    
	}

	close $seq_fh;

	# Test empty file
	if (-z $seqfile) {
	    print STDERR "Error, empty sequence file...\n";
	}

	# print STDERR "Running GeneID, with this command:\n";
	# print STDERR "$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile \n";

        my $geneid_output = qx/$_geneid_dir\/$_geneid_bin $_geneid_args $seqfile/;
        
	# Comment this thine if you want to keep the file...
	unlink $seqfile;

        if (defined $geneid_output) {
		return $geneid_output;
	}	
	else {
		# What else better to return ??
		return undef;
	}
}


=head2 SGP2_call

 Title   : SGP2_call
 Usage   : $report = SGP2_call (@params);
         : 
         : ## where @params are,
         : @params = ('arg1'  => "WKRPPEICENPRFIIGGANRTDIAAIACLTLNERL",
         :            'arg2'  => "query 1", ## optional
         :            'arg3'  => "nr");     ## optional (default: nr)
 Returns : Devuelve un string que contiene el resultado de la ejecuci�n.

=cut

sub SGP2_call {
        my %args = @_; 
        # relleno los parametros por defecto SGP2_call (nucleotide  => $nucleotide, seqIdentifier => $sequenceIdentifier);
        my $sequences          = $args{sequences} || undef;
        my $tblastx_output     = $args{tblastx_output} || undef;
	my $format             = $args{format} || "";
	my $parameters         = $args{parameters} || undef;

	# Get the parameters
	# No one yet!

        # Llama a SGP2 localmente
        my $_sgp2_dir  = "/home/ug/gmaster/sgp2/sgp2_2003/";
	$_sgp2_dir  = $ENV{SGP2};
        my $_sgp2_bin  = "bin/sgp2";
        my $_sgp2_args = "";
	
        if ($format eq "GFF") {
	    $_sgp2_args .= "-g G";
        }
	
	# Generate a temporary file locally with the sequence in FASTA format
	# locally, ie not on a NFS mounted directory, for speed sake

	# my ($fh, $seqfile) = tempfile("/tmp/moby/SGP2_Sequence.XXXXXX", UNLINK => 1);
	# close($fh);

	my @seqIds = keys (%$sequences);
	foreach my $sequenceIdentifier (@seqIds) {

	    my $nucleotide = $sequences->{$sequenceIdentifier};

	    # bioperl object
	    
	    my $seqobj = Bio::Seq->new (
					-display_id => $sequenceIdentifier,
					-sequence   => $nucleotide
					);
	    
	    # Bioperl sequence factory
	    
	    # my $sout = Bio::SeqIO->new (
	    # 			    -file   => ">$seqfile",
	    #			    -format => 'fasta'
	    #			    );
	    # $sout->write_seq ($seqobj);
	    
	}

	my $seqfile = "/home/ug/arnau/projects/sgp2/sgp2/samples/Hsap_BTK.msk.fa";

	# TBLASTX Output File

	# Generate a temporary file locally with the TBLASTX Output
	# locally, ie not on a NFS mounted directory, for speed sake

	# my ($fh, $tblastx_output_file) = tempfile("/tmp/moby/SGP2_TBLASTX.XXXXXX", UNLINK => 1);
	# close($fh);

	# qx/echo "$tblastx_output" > $tblastx_output_file/;

	my $tblastx_output_file = "/home/ug/arnau/projects/sgp2/sgp2/samples/Hsap_BTK.tbx";

	# Test empty files

	if (-z $seqfile) {
	    print STDERR "Error, empty sequence file...\n";
	}

	if (-z $tblastx_output_file) {
	    print STDERR "Error, empty tblastx output file...\n";
	}

	# print STDERR "Running SGP2, with this command:\n";
	# print STDERR "$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file\n";

        my $sgp2_output = qx/$_sgp2_dir\/$_sgp2_bin $_sgp2_args -1 $seqfile -t $tblastx_output_file/;
        
        if (defined $sgp2_output) {
		return $sgp2_output;
	}	
	else {
		# What else better to return ??
		return undef;
	}
}

1;

__END__
