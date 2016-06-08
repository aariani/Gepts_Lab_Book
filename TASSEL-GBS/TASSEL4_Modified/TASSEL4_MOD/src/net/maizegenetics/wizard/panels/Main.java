//package net.maizegenetics.wizard.panels;
//
//import java.awt.event.MouseListener;
//import net.maizegenetics.wizard.WizardPanelDescriptor;
//import net.maizegenetics.wizard.Wizard;
////import net.maizegenetics.wizard.*;
//
//import javax.swing.*;
//
//public class Main {
//
//    public static void main(String[] args) {
//
//        Wizard wizard = new Wizard();
//        wizard.getDialog().setTitle("TASSEL Wizard (In development)");
//
//        WizardPanelDescriptor descriptor1 = new TestPanel1Descriptor();
//        wizard.registerWizardPanel(TestPanel1Descriptor.IDENTIFIER, descriptor1);
//
//        WizardPanelDescriptor descriptor2 = new TestPanel2Descriptor();
//        wizard.registerWizardPanel(TestPanel2Descriptor.IDENTIFIER, descriptor2);
//
//        WizardPanelDescriptor descriptor3 = new TestPanel3Descriptor();
//        wizard.registerWizardPanel(TestPanel3Descriptor.IDENTIFIER, descriptor3);
//
//        wizard.setCurrentPanel(TestPanel1Descriptor.IDENTIFIER);
//
//        int ret = wizard.showModalDialog();
//
//        System.out.println("Dialog return code is (0=Finish,1=Cancel,2=Error): " + ret);
//        System.out.println("Second panel selection is: " +
//            (((TestPanel2)descriptor2.getPanelComponent()).getRadioButtonSelected()));
//
//        //System.exit(0);
//
//    }
//
//}
