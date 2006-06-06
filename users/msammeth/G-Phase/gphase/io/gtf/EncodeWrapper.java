package gphase.io.gtf;

import gphase.algo.ASAnalyzer;
import gphase.model.ASMultiVariation;
import gphase.model.ASVariation;
import gphase.model.AbstractRegion;
import gphase.model.DefaultRegion;
import gphase.model.DirectedRegion;
import gphase.model.Exon;
import gphase.model.Gene;
import gphase.model.Graph;
import gphase.model.Species;
import gphase.model.Transcript;
import gphase.tools.Arrays;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.Collection;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Vector;

import com.p6spy.engine.logging.appender.StdoutLogger;

public class EncodeWrapper extends GTFWrapper {
	
	public Graph getGraph(boolean encode) {
		try {
			read();
		} catch (Exception e) {
			e.printStackTrace(); 
		}
		return assemble(encode);		// <===== check ENCODE here !!!
	}
	
	public static void main(String[] args) {
		//"encode/44regions_genes_CHR_coord.gtf"
		//"encode/EnsemblGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_fromUCSC_inENCODEonly.gtf"
		// "encode/RefSeqGenes_fromUCSC.inENCODE.gtf"
		// "encode/RefSeqGenes_fromUCSC.gtf"
		// "encode/EnsemblGenes_all_fromENSEMBL.gtf"
		// "encode/RefseqGenes_all_fromENSEMBL.gtf"
		String fName= "encode/RefseqGenes_all_fromENSEMBL.gtf";
		EncodeWrapper myWrapper= new EncodeWrapper(new File(fName).getAbsolutePath()); // testGTF.gtf
		try {
			myWrapper.read();
		} catch (Exception e) {
			e.printStackTrace(); 
		}
		boolean encode= false;
		if (fName.equals("encode/44regions_genes_CHR_coord.gtf"))
			encode= true;
		
		long t0= System.currentTimeMillis();
		Graph g= myWrapper.assemble(encode);		// <===== check ENCODE here !!!
		System.out.println(g.getTranscriptCount());
		g.filterNonCodingTranscripts();
		System.out.println(g.getTranscriptCount());

		
		//ASAnalyzer.outputASStatistics(g);
		ASAnalyzer.analyze1_transcripts_clusters(g);
		
//		ASVariation[][] as= ASAnalyzer.determineVariations(g, (PrintStream) null);
//		for (int i = 0; i < as.length; i++) {
//			System.out.println(as[i][0].toString()+"\t"+as[i][0].toBitString());
//		}
		
//		ASAnalyzer.lengthVariationModulo(g);
//		long t1= System.currentTimeMillis();
//		System.out.println("build up"+(t1- t0));
//		g.getASVariations(ASMultiVariation.FILTER_NONE);
//		long t2= System.currentTimeMillis();
//		System.out.println("filter: "+ (t2- t1));
//		g.getASVariations(ASMultiVariation.FILTER_REDUNDANT);
//		long t3= System.currentTimeMillis();
//		System.out.println("filter redundant: "+(t3- t2));
		
			// analyze1
//		System.out.println("clusters: "+g.getGenes().length+" transcripts: "+g.getTranscriptCount());
//		ASAnalyzer.filterSingleExonClusters(g);
//		System.out.println("clusters: "+g.getGenes().length+" transcripts: "+g.getTranscriptCount());
//		System.out.println("clusters w as: "+g.getClustersWithAS());
		
		//g.filterNonCodingTranscripts();
		//System.out.println(ASAnalyzer.getPercASSpliceSites(g));
		//ASVariation[][] classes= ASAnalyzer.determineVariations(g, "isTrue");
		//ASAnalyzer.lengthVariation(g);
		//ASAnalyzer.getSpecificVariations(g, System.out);
		//ASAnalyzer.getVincenco(g);
		//ASAnalyzer.getAllASVariation(g);
		//ASAnalyzer.determineVariations(g, System.out);
//		System.out.println(g.getTranscriptCount());
		
		
//		ASAnalyzer.analyzeGene(g, "RNPC2", "isTrue", ASMultiVariation.FILTER_REDUNDANT);
//		int cnt= 0, not= 0;
//		BufferedWriter buffy= new BufferedWriter(new FileWriter)
//		for (int i = 0; i < classes.length; i++) {
//			for (int j = 0; j < classes[i].length; j++) {
//				classes[i][j].outputDetail(System.out);
//			}
//		}
//		System.out.println(cnt+","+not+" ("+(cnt+not)+")");
		
		
//		ASVariation[][] classes= null;
//		ASVariation[] events= ASAnalyzer.getVariation("(1=4^ // 2=3^)", classes);
//		int cnt= 0;
//		for (int i = 0; i < classes.length; i++) {
//			for (int j = 0; events!= null&& j < events.length; j++) {
//				if (events[j].includesFirstSpliceSite()) {
//					events[j].outputDetail(System.out);
//				}
			//}
//		}
//		System.out.println("--");
				
		
//			events= ASAnalyzer.getVariation("(1=3^ // 2=4^)", classes);
////			int cnt= 0;
////			for (int i = 0; i < classes.length; i++) {
//				for (int j = 0; j < events.length; j++) {
////					if (events[j].includesFirstSpliceSite()) {
//						events[j].outputDetail(System.out);
//					}
				//}
//			}
//			System.out.println(cnt);
		//ASVariation[][] classes2= ASAnalyzer.determineVariations(g, "isNotAtAllCoding", true);
//		ASAnalyzer.diff(classes,classes2);
//		classes= ASAnalyzer.getVariation("(1^ // 2^)", classes);
//		classes= ASAnalyzer.getVariation("( // 1^2=)", classes);
//		for (int i = 0; i < classes.length; i++) {
//			for (int j = 0; j < 10; j++) {
//				classes[i][j].outputDetail(System.out);
//			}
//		}
		
		//ASAnalyzer.determineVariations(g, "isPartiallyCoding");
		
		//ASAnalyzer.analyzeGene(g, "")
		//ASAnalyzer.lowRepresented(g, 10);
		//ASAnalyzer.debug(g, 10);
		//ASAnalyzer.lengthPlot(g);
		//ASAnalyzer.getVariation("( // 1=2^3=4^5=6^7=8^9=10^11=12^13=14^15=16^17=18^19=20^21=22^23=24^25=26^27=28^)",
		//		g);

	}
	
	public EncodeWrapper(String absFName) {
		super(absFName);
	}
	
	GTFObject createGTFObject(){
		return new France();
	}

	public France[] getFranceObj() {
		return (France[]) getGtfObj();
	}
	
	Vector getClusters(HashMap transHash) {
		Vector trans= new Vector(transHash.values());
		HashMap clusters= new HashMap();	// transcriptID x (vector x vector)
		for (int i = 0; i < trans.size(); i++) {
			for (int j = i+ 1; j < trans.size(); j++) {
				if (overlap((Vector) trans.elementAt(i), (Vector) trans.elementAt(j))) {
					Vector v1= (Vector) trans.elementAt(i);
					Vector v2= (Vector) trans.elementAt(j);
					Vector vc1= (Vector) clusters.remove(((France) v1.elementAt(0)).getTranscriptID());
					// remove from all other transcript ids in the cluster in hash
					if (vc1== null) {
						vc1= new Vector();
						vc1.add(v1);
					}
					Vector vc2= (Vector) clusters.remove(((France) v2.elementAt(0)).getTranscriptID());
					if (vc2== null) {
						vc2= new Vector();
						vc2.add(v2);
					}
					Vector result= (Vector) Arrays.merge(vc1, vc2);	// doppelte eintraege moeglich??
					for (int k = 0; k < result.size(); k++) 
						clusters.put(((France) ((Vector) result.elementAt(k)).elementAt(0)).getTranscriptID(),
								result);
				}
					
			}
		}
		return trans;
	}

	boolean overlap (Vector trans1, Vector trans2) {
//		for (int i = 0; i < trans1.size(); i++) {
//			France f1= (France) trans1.elementAt(i);
//			if (!f1.isExon())
//				continue;
//			for (int j = 0; j < trans2.size(); j++) {
//				France f2= (France) trans2.elementAt(j);
//				if (!f2.isExon())
//					continue;
//				DefaultDirectedRegion d1= new DefaultDirectedRegion(f1.getStrand(), f1.getStart(), f1.getEnd());
//				DefaultDirectedRegion d2= new DefaultDirectedRegion(f1.getStrand(), f2.getStart(), f2.getEnd());
//				if (d1.overlaps(d2))
//					return true;
//			}
//		}
		return false;
	}
	
	HashMap getChromosomes(HashMap transGTF) {
		Iterator iter= transGTF.values().iterator();
		HashMap chrHash= new HashMap();
		while (iter.hasNext()) {
			Vector gtfsVec= (Vector) iter.next();
			France o= (France) (gtfsVec).elementAt(0);
			HashMap tHash= (HashMap) chrHash.remove(o.getChromosome());
			if (tHash== null)
				tHash= new HashMap();
			Vector v= (Vector) tHash.remove(o.getTranscriptID());
			if (v== null)
				v= new Vector();
			for (int i = 0; i < gtfsVec.size(); i++) 
				v.add(gtfsVec.elementAt(i));
			tHash.put(o.getTranscriptID(), v);
			chrHash.put(o.getChromosome(), tHash);
		}
		return chrHash;
	}
	
	Graph assemble(boolean encode) {
		
		Species spec= new Species("human");

			// cluster
		HashMap hash= getGroups("transcript_id", getGtfObj());	// cluster for genes?
		HashMap chrHash= getChromosomes(hash);
		
			// construct transcripts
		Collection co= ((Collection) chrHash.keySet());
		String[] keys= new String[co.size()];
		Iterator iter= co.iterator();
		int x= 0;
		while(iter.hasNext()) 
			keys[x++]= (String) iter.next();
		
		HashMap chr2Hash= new HashMap(chrHash.size());
		for (int i = 0; i < keys.length; i++) {	// chromosomes
			String chrID= keys[i];
			HashMap tHash= (HashMap) chrHash.get(chrID);
			Collection co2= ((Collection) tHash.keySet());
			String[] tkeys= new String[co2.size()];
			Iterator iter2= co2.iterator();
			x= 0;
			while (iter2.hasNext())					
				tkeys[x++]= (String) iter2.next();
			HashMap t2Hash= new HashMap(tHash.size());
			chr2Hash.put(chrID, t2Hash);
			for (int j = 0; j < tkeys.length; j++) {	// transcripts
				String tID= tkeys[j];
				GTFObject[] gtfs= (GTFObject[]) Arrays.toField(tHash.get(tID));	// gtf entries for 1 transcript
				France ff= (France) gtfs[0];
				if (encode&& !ff.getSource().contains("VEGA"))
					continue;
				Transcript transcript= new Transcript(tID);
				transcript.setStrand(ff.getStrand());
				for (int k = 0; k < gtfs.length; k++) {		// exons 
					France f= (France) gtfs[k];
					if (f.isExon()) 
						transcript.setBoundaries(new Exon(transcript, f.getExonID(), f.getStart(), f.getEnd()));
				}
				t2Hash.put(tID, transcript);	// fill tHash with transcripts
			}
			
		}
		
			// cluster
		HashMap gHash= new HashMap();
		Comparator compi= new AbstractRegion.PositionComparator();
		for (int i = 0; i < keys.length; i++) {	// chromosomes
			String chrID= keys[i];
			HashMap t2Hash= (HashMap) chr2Hash.get(chrID);
			Object[] transcripts= t2Hash.values().toArray();
			java.util.Arrays.sort(transcripts, compi);
			Transcript[] t= new Transcript[transcripts.length];
			for (int j = 0; j < t.length; j++) 
				t[j]= (Transcript) transcripts[j];
			Transcript[][] loci= clusterTranscripts(t);
			for (int j = 0; j < loci.length; j++) {
				String gID= Gene.getUniqueID();
				Gene locus= new Gene(spec, gID);
				locus.setStrand(loci[j][0].getStrand());
				locus.setChromosome(chrID);
				for (int k = 0; k < loci[j].length; k++) { // transcripts
					loci[j][k].setGene(locus);
					Vector v= (Vector) ((HashMap) chrHash.get(chrID)).get(loci[j][k].getTranscriptID());
					for (int m = 0; m < v.size(); m++) {
						France f= (France) v.elementAt(m);
						if (f.isExon())
							loci[j][k].addExon(new Exon(loci[j][k], f.getExonID(), f.getStart(), f.getEnd()));
						else if (f.isCDS())
							loci[j][k].addCDS(f.getStart(), f.getEnd());
					}
					locus.addTranscript(loci[j][k]);
				}
				gHash.put(gID, locus);
			}
		}		
		
			// build graph
		iter= gHash.values().iterator();
		Graph g= new Graph();
		g.addSpecies(spec);
		while (iter.hasNext()) 
			g.addGene((Gene) iter.next());
		return g;
	}
	
	

	private Transcript[][] clusterTranscripts(AbstractRegion[] regions) {
		
		int max= Integer.MIN_VALUE;
		Vector clusters= new Vector();
		Vector v= null;
		for (int i = 0; i < regions.length; i++) {

			if (regions[i].getStart()> max) {
				if (v!= null)
					clusters.add(v);
				v= new Vector();
			} 
			v.add(regions[i]);
			if (regions[i].getEnd()> max)
				max= regions[i].getEnd();
		}
		if (v!= null)
			clusters.add(v);
		
		return (Transcript[][]) Arrays.toField(clusters);
	}

	Graph assemble_stable(boolean encode) {
		
		Species spec= new Species("human");
			// cluster
		HashMap gHash= getGroups("gene_id", getGtfObj());	// cluster for genes
		
			// infer objects
		Collection co= ((Collection) gHash.keySet());
		String[] keys= new String[co.size()];
		Iterator iter= co.iterator();
		int x= 0;
		while(iter.hasNext()) 
			keys[x++]= (String) iter.next();
		
		for (int i = 0; i < keys.length; i++) {	// genes
			String key= keys[i];
			Gene gene= new Gene(spec, key);
			GTFObject[] gtfs= (GTFObject[]) Arrays.toField(gHash.remove(key));	// transcripts in gene
			France ff= (France) gtfs[0];
			gene.setStrand(ff.getStrand());
			gene.setChromosome(ff.getChromosome());
			gene.setConfidence(ff.getSource());
			HashMap tHash= getGroups("transcript_id", gtfs);
			Collection co2= ((Collection) tHash.keySet());
			String[] tkeys= new String[co2.size()];
			Iterator iter2= co2.iterator();
			x= 0;
			while (iter2.hasNext())
				tkeys[x++]= (String) iter2.next(); 
			for (int j = 0; j < tkeys.length; j++) {	// transcripts
				String tID= tkeys[j];
				Transcript transcript= new Transcript(gene, tID);
				Vector v= (Vector) tHash.remove(tID);
				for (int k = 0; k < v.size(); k++) {		// exons 
					France f= (France) v.elementAt(k);
					if (f.isExon()) 
						transcript.addExon(new Exon(transcript, f.getExonID(), f.getStart(), f.getEnd()));
					else if (f.isCDS())
						transcript.addCDS(f.getStart(), f.getEnd());
				}
				gene.addTranscript(transcript);
			}
			if (!encode|| gene.getConfidence().contains("VEGA"))
				gHash.put(key, gene);	// re-fill hash
		}
		
			// build graph
		iter= gHash.values().iterator();
		Graph g= new Graph();
		g.addSpecies(spec);
		while (iter.hasNext()) 
			g.addGene((Gene) iter.next());
		return g;
	}
	
	
	void assemble_old() {
		
//		Species spec= new Species("human");
//		Vector clusters= getClusters(transHash);	// Vector x Vector x Vector
//		
//		Collection c= transHash.values();
//		HashMap tHash= new HashMap(c.size());
//		HashMap gHash= new HashMap(c.size()/ 2);
//		Iterator it= clusters.iterator();
//		while (it.hasNext()) {	// iterate all clusters
//			Vector v= (Vector) it.next();
//			Cluster currCluster= new Cluster(((France) ((Vector) v.elementAt(0)).elementAt(0)).getChromosome(), spec);
//			for (int i = 0; i < v.size(); i++) {
//				Vector vTrans= (Vector) v.elementAt(i);
//				for (int j = 0; j < vTrans.size(); j++) {
//					France f= (France) vTrans.elementAt(j);	// exon, cds, ...
//					Transcript t= (Transcript) tHash.get(f.getTranscriptID());	// get corresponding transcript
//					if (t== null) {
//						t= new Transcript(f.getTranscriptID(),f.getStrand(),currCluster);
//						t.setSource(f.getSource());
//						Gene g= (Gene) gHash.get(f.getGeneID());
//						if (g== null) {
//							g= new Gene(f.getGeneID());
//							g.setAlias(f.getGeneAlias());
//							gHash.put(f.getGeneID(), g);
//						}
//						t.addGene(g);
//						tHash.put(f.getTranscriptID(), t);
//					}
////					if (f.isCDS())
////						t.addCDS(f.getExonID(), f.getStart(), f.getEnd());
//				}
//			}
//		}
	}
	
	private HashMap getGroups(String id, GTFObject[] obj) {
		
		HashMap hash= new HashMap();
		for (int i = 0; i < obj.length; i++) {
			Vector tAttrib= (Vector) hash.get(obj[i].getAttribute(id));
			if (tAttrib== null) {
				tAttrib= new Vector();
				hash.put(obj[i].getAttribute(id), tAttrib);
			}
			tAttrib.add(obj[i]);
		}
		
		return hash;
	}

	private HashMap getTranscriptInfo() {
		
		HashMap tHash= new HashMap();
		France[] france= getFranceObj();
		for (int i = 0; i < france.length; i++) {
			if (!france[i].isExon()&& !france[i].isCDS())	// collect exons and cds
				continue;
			Vector tAttrib= (Vector) tHash.get(france[i].getTranscriptID());
			if (tAttrib== null) {
				tAttrib= new Vector();
				tHash.put(france[i].getTranscriptID(), tAttrib);
			}
			tAttrib.add(france[i]);
		}
		
		return tHash;
	}
	
}
