/*
 Copyright (C) 2003 EBI, GRL

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

package org.ensembl.util;

import java.util.HashSet;
import java.util.Iterator;

/**
 * Set for storing unique longs.
 * 
 * @author <a href="mailto:craig@ebi.ac.uk">Craig Melsopp </a>
 */
public class LongSet extends HashSet {

	private static final long serialVersionUID = 1L;

	protected long[] result = null;

	public void addAll(long[] ls) {
		for (int i = 0; i < ls.length; i++) 
			add(ls[i]);
	}
	
	public void add(long l) {
		add(new Long(l));
		result = null;
	}

	public long[] to_longArray() {
		if (result == null) {
			final int size = size();
			result = new long[size];
			Iterator iter = iterator();
			for (int i = 0; i < size; i++)
				result[i] = ((Long) iter.next()).longValue();
		}
		return result;
	}

}
