// Decompiled by DJ v3.9.9.91 Copyright 2005 Atanas Neshkov  Date: 12/05/2006 07:41:05 PM
// Home Page : http://members.fortunecity.com/neshkov/dj.html  - Check often for new version!
// Decompiler options: packimports(3) fieldsfirst nonlb 
// Source File Name:   GenePane.java

import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import javax.swing.JColorChooser;
import javax.swing.JFrame;

private final class GenePane$2
    implements ActionListener {

    public void actionPerformed(ActionEvent e) {
        java.awt.Color newColor = JColorChooser.showDialog(null, "Choose UTR Color", GenePane.UTRCOLOR1);
        if(newColor == null) {
            return;
        } else {
            GenePane.UTRCOLOR1 = newColor;
            optionPane.setVisible(false);
            tion.access._mth0(tion.this).repaint();
            return;
        }
    }

    GenePane$2() {
    }
}