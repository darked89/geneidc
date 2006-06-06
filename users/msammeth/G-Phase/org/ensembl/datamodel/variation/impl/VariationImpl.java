/*
    Copyright (C) 2001 EBI, GRL

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; either
    version 2.1 of the License, or (at your option) any later version.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
 */
package org.ensembl.datamodel.variation.impl;

import java.util.ArrayList;
import java.util.List;

import org.ensembl.datamodel.impl.PersistentImpl;
import org.ensembl.datamodel.variation.Allele;
import org.ensembl.datamodel.variation.ValidationState;
import org.ensembl.datamodel.variation.Variation;
import org.ensembl.driver.variation.VariationDriver;
import org.ensembl.util.StringUtil;

/**
 * A sequence variation.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public class VariationImpl extends PersistentImpl implements Variation {

  /**
   * Used by the (de)serialization system to determine if the data 
   * in a serialized instance is compatible with this class.
   *
   * It's presence allows for compatible serialized objects to be loaded when
   * the class is compatible with the serialized instance, even if:
   *
   * <ul>
   * <li> the compiler used to compile the "serializing" version of the class
   * differs from the one used to compile the "deserialising" version of the
   * class.</li>
   *
   * <li> the methods of the class changes but the attributes remain the same.</li>
   * </ul>
   *
   * Maintainers must change this value if and only if the new version of
   * this class is not compatible with old versions. e.g. attributes
   * change. See Sun docs for <a
   * href="http://java.sun.com/j2se/1.4.2/docs/guide/serialization/">
   * details. </a>
   *
   */
  private static final long serialVersionUID = 1L;



  private List validationStates = new ArrayList();


  private List synonyms = new ArrayList();


  private List alleles = new ArrayList();


  private VariationDriver vdriver;


  private String name;


  /**
   * This empty constructor is designed to be used by 
   * VariationFeatureImpl only.
   */
  VariationImpl() {
    
  }
  
  
  public VariationImpl(VariationDriver vdriver) {
    this.vdriver = vdriver;
  }



  /**
   * @see org.ensembl.datamodel.variation.Variation#getSynonyms()
   */
  public List getSynonyms() {
    return synonyms;
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#getSynonymSources()
   */
  public List getSynonymSources() {
    // TODO Auto-generated method stub
    return null;
  }

  /**
   * @see org.ensembl.datamodel.variation.Variation#addSynonym(java.lang.String)
   */
  public void addSynonym(String synonym) {
    synonyms.add(synonym);
    
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#getValidationStates()
   */
  public List getValidationStates() {
    return validationStates;
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#addValidationState(org.ensembl.datamodel.variation.ValidationState)
   */
  public void addValidationState(ValidationState state) {
    validationStates.add(state);
    
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#getSource()
   */
  public String getSource() {
    // TODO Auto-generated method stub
    return null;
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#getAlleles()
   */
  public List getAlleles() {
    return alleles;
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#addAllele(org.ensembl.datamodel.variation.Allele)
   */
  public void addAllele(Allele allele) {
    alleles.add(allele);
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#getFivePrimeFlankingSeq()
   */
  public String getFivePrimeFlankingSeq() {
    // TODO Auto-generated method stub
    return null;
  }

  /* (non-Javadoc)
   * @see org.ensembl.datamodel.variation.Variation#getThreePrimeFlankingSeq()
   */
  public String getThreePrimeFlankingSeq() {
    // TODO Auto-generated method stub
    return null;
  }

  /**
   * @see org.ensembl.datamodel.variation.Variation#getName()
   */
  public String getName() {
    return name;
  }

  /**
   * @see org.ensembl.datamodel.variation.Variation#setName()
   */
  public void setName(String name) {
    this.name=name;    
  }

  
  /**
   * @see org.ensembl.datamodel.impl.PersistentImpl#toString()
   */
  public String toString() {
    StringBuffer buf = new StringBuffer("[");
    buf.append(super.toString());
    buf.append("vdiver=").append(StringUtil.setOrUnset(vdriver));
    buf.append(", name=").append(name);
    buf.append(", synonyms=").append(StringUtil.sizeOrUnset(synonyms));
    buf.append(", validationStates=").append(StringUtil.sizeOrUnset(validationStates));
    return buf.append("]").toString();
  }
}
