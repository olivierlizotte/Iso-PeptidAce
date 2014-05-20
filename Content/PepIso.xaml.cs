using System;
using System.Collections.Generic;
using System.Linq;
using System.Text;
using System.Threading.Tasks;
using System.Windows;
using System.Windows.Controls;
using System.Windows.Data;
using System.Windows.Documents;
using System.Windows.Input;
using System.Windows.Media;
using System.Windows.Media.Imaging;
using System.Windows.Navigation;
using System.Windows.Shapes;
using PeptidAce.Iso.Methods;
using PeptidAce.Utilities.Interfaces;

namespace PeptidAce.ModernUI.Content
{
    /// <summary>
    /// Interaction logic for PepIso.xaml
    /// </summary>
    public partial class PepIso : UserControl
    {
        public static PositionnalIsomerSolver solver;
        public static double precursorMassTolPPm = 8;
        public static double productMassTolPPm   = 0.05;
        public static int nbMinFragments = 5;
        public static int nbMaxFragments = 5;
        public string[] peptideFiles;
        public string[] mixedFiles;
        public string   fastaFile;

        public PepIso()
        {
            InitializeComponent();            
        }

        private void Button_Click_OpenFilePeptide(object sender, RoutedEventArgs e)
        {            
            bool? dResult = false;
            Microsoft.Win32.OpenFileDialog open_dlg = new Microsoft.Win32.OpenFileDialog();
            open_dlg.Title = "Select Raw files describing individual peptides";
            //open_dlg.Filter = "*";
            open_dlg.Multiselect = true; 
            dResult = open_dlg.ShowDialog();
            if(open_dlg.FileNames != null)
            {
                string text = "";
                peptideFiles = open_dlg.FileNames;
                foreach(string file in open_dlg.FileNames)
                    text += file + "\n";
                PeptideRawFiles.Text = text;
                PeptideRawFiles.Foreground = Brushes.Black;
            }
        }

        private void Button_Click_OpenFileMixed(object sender, RoutedEventArgs e)
        {
            bool? dResult = false;
            Microsoft.Win32.OpenFileDialog open_dlg = new Microsoft.Win32.OpenFileDialog();
            open_dlg.Title = "Select Raw files of mixed spectrum";
            //open_dlg.Filter = "*";
            open_dlg.Multiselect = true;
            dResult = open_dlg.ShowDialog();
            if (open_dlg.FileNames != null)
            {
                string text = "";
                mixedFiles = open_dlg.FileNames;
                foreach (string file in open_dlg.FileNames)
                    text += file + "\n";
                MixedRawFiles.Text = text;
                MixedRawFiles.Foreground = Brushes.Black;
            }
        }

        private void Button_Click(object sender, RoutedEventArgs e)
        {

            bool? dResult = false;
            Microsoft.Win32.OpenFileDialog open_dlg = new Microsoft.Win32.OpenFileDialog();
            open_dlg.Title = "Select the fasta file with your sequences of interest";
            //open_dlg.Filter = "*";
            open_dlg.Multiselect = false;
            dResult = open_dlg.ShowDialog();
            if (open_dlg.FileNames != null)
            {
                fastaFile = open_dlg.FileName;
                FastaFile.Text = fastaFile;
                FastaFile.Foreground = Brushes.Black;
            }
        }

        public Task runningTask = null;
        public ConSolBasic conSol = null;

        public void AddTextOutput(string msg)
        {
            TextConsol.Dispatcher.Invoke(
                System.Windows.Threading.DispatcherPriority.Normal,
                new Action(
                    delegate()
                    {
                        TextConsol.Text += "\n" + msg;
                        ScrollText.ScrollToBottom();
                    }
                ));
        }

        private bool CheckParams()
        {
            string msg = "";
            if(precursorMassTolPPm <= 0)
                msg += "\nPrecursor Mass Tolerance needs to be above 0. Check Advanced Settings";

            if (productMassTolPPm <= 0)
                msg += "\nFragment Mass Tolerance needs to be above 0. Check Advanced Settings";

            
            if (nbMinFragments <= 1)
                msg += "\nThe Minimum number of Fragments has to be bigger than 1. Check Advanced Settings";

            if (nbMaxFragments > 100)
                msg += "\nThe Maximum number of Fragments has to be smaller than 100. Check Advanced Settings";

            if(nbMinFragments > nbMaxFragments)
                msg += "\nThe minimum number of fragments cannot be smaller than the maximum number of fragments. Check Advanced Settings";

            if (peptideFiles == null || peptideFiles.Length < 2)
                msg += "\nAt least two RAW files with individual peptides are necessary";

            if (mixedFiles == null || mixedFiles.Length < 1)
                msg += "\nAt least one mixed spectrum RAW file is necessary";

            if (string.IsNullOrEmpty(fastaFile))
                msg += "\nPlease select a Fasta file with the sequence you wich to search for. This can be the list of peptide sequences or your protein of interest. Keep this file small for faster searches.";
            
            if (!string.IsNullOrEmpty(msg))
            {
                AddTextOutput(msg);
                return false;
            }
            else return true;
        }
        private void Button_Click_1(object sender, RoutedEventArgs e)
        {
            if (conSol == null)
                conSol = new ConSolBasic(AddTextOutput);


            //PeptidAce.Iso.UnitTests.StatsMaker.ProjectMerge(conSol);
            //PeptidAce.Iso.UnitTests.PisTest.Uptimize();

            if (ButtonRun.IsEnabled && CheckParams())
            {
                ButtonRun.IsEnabled = false;
                TextConsol.Text = "Launching Peptide Isomer solver ... ";

                runningTask = new Task(() =>
                {
                    PositionnalIsomerSolver newSolver = new PositionnalIsomerSolver();
                    try
                    {
                        newSolver.precTolPpm = precursorMassTolPPm;
                        newSolver.prodTolDa = productMassTolPPm;
                        newSolver.nbMinFragments = nbMinFragments;
                        newSolver.nbMaxFragments = nbMaxFragments;
                        newSolver.Solve(peptideFiles, mixedFiles, fastaFile, Utilities.vsCSV.GetFolder(mixedFiles[0]), conSol);

                        solver = newSolver;
                    }
                    catch(Exception ex)
                    {
                        AddTextOutput(ex.Message);
                        AddTextOutput(ex.StackTrace);
                    }

                    try
                    {
                        ButtonRun.Dispatcher.Invoke(
                            System.Windows.Threading.DispatcherPriority.Normal,
                            new Action(
                                delegate()
                                {
                                    ButtonRun.IsEnabled = true;
                                    AddTextOutput("\nDone! Your results are available in 'csv' files in these folders :\n" +
                                                  "\n" + Utilities.vsCSV.GetFolder(mixedFiles[0]) + "Identifications" +
                                                  "\n" + Utilities.vsCSV.GetFolder(mixedFiles[0]) + "Combined" +
                                                  "\n" + Utilities.vsCSV.GetFolder(mixedFiles[0]) + "Individual");

                                    ModernUI.Pages.Home.UpdateLinks(newSolver.SpikedSamples, newSolver.mixedPrecursors);
                                }
                            ));
                    }
                    catch (Exception ex)
                    {
                        AddTextOutput(ex.Message);
                        AddTextOutput(ex.StackTrace);
                    }
                });
                runningTask.Start();
            }
        }

        private void PeptideRawFiles_Drop(object sender, DragEventArgs e)
        {
                // If the DataObject contains string data, extract it. 
                if (e.Data.GetDataPresent(DataFormats.FileDrop, true))
                {

                    peptideFiles = e.Data.GetData(DataFormats.FileDrop, false) as string[];
                    PeptideRawFiles.Text = "";
                    if(peptideFiles.Length > 0)
                    {
                        PeptideRawFiles.Text = peptideFiles[0];
                        for (int i = 1; i < peptideFiles.Length; i++)
                            PeptideRawFiles.Text += "\n" + peptideFiles[i];
                    }
                }
        }

        private void PeptideRawFiles_DragEnter(object sender, DragEventArgs e)
        {
            e.Handled = true;
            e.Effects = DragDropEffects.Copy;
        }

        private void StackPanel_DragEnter(object sender, DragEventArgs e)
        {
            e.Handled = true;
            e.Effects = DragDropEffects.Copy;
        }

        private void StackPanel_Drop(object sender, DragEventArgs e)
        {
            // If the DataObject contains string data, extract it. 
            if (e.Data.GetDataPresent(DataFormats.FileDrop, true))
            {

                mixedFiles = e.Data.GetData(DataFormats.FileDrop, false) as string[];
                MixedRawFiles.Text = "";
                if (mixedFiles.Length > 0)
                {
                    MixedRawFiles.Text = mixedFiles[0];
                    for (int i = 1; i < mixedFiles.Length; i++)
                        MixedRawFiles.Text += "\n" + mixedFiles[i];
                }
            }
        }

        private void StackPanel_DragEnter_1(object sender, DragEventArgs e)
        {
            e.Handled = true;
            e.Effects = DragDropEffects.Copy;
        }

        private void StackPanel_Drop_1(object sender, DragEventArgs e)
        {
            // If the DataObject contains string data, extract it. 
            if (e.Data.GetDataPresent(DataFormats.FileDrop, true))
            {
                string[] files = e.Data.GetData(DataFormats.FileDrop, false) as string[];
                if (files.Length == 1)
                {
                    fastaFile = files[0];
                    FastaFile.Text = fastaFile;
                }
            }
        }
    }
}
