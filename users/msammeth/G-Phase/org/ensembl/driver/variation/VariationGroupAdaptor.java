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
package org.ensembl.driver.variation;

import java.util.List;

import org.ensembl.datamodel.variation.Variation;
import org.ensembl.datamodel.variation.VariationGroup;
import org.ensembl.driver.Adaptor;
import org.ensembl.driver.AdaptorException;

/**
 * Retrieves Variation groups from an ensembl database.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp</a>
 *
 */
public interface VariationGroupAdaptor extends Adaptor {

  final static String TYPE = "variation_group";

  /**
   * Fetch Variation group with specified internal ID.
   * @return specified variation group or null if none found.
   */
  VariationGroup fetch(long internalID) throws AdaptorException;

  /**
   * Fetch Variation group with specified name.
   * @return specified variation group or null if none found.
   */
  VariationGroup fetch(String name) throws AdaptorException;

  /**
   * Fetch all variation groups that contain the variation.
   * @return zero or more VariationGroups.
   */
  List fetch(Variation variation) throws AdaptorException;
}
