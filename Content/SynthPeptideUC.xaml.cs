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
using System.ComponentModel;
using System.Collections.ObjectModel;

namespace PeptidAce.ModernUI.Content
{
    /// <summary>
    /// Interaction logic for RatioUC.xaml
    /// </summary>
    public partial class SynthPeptideUC : UserControl
    {
        ObservableCollection<PeptideDataItem> data = new ObservableCollection<PeptideDataItem>()
            {
                //new PieRatioDataItem() { Title = "s1", Value = 10 },
                //new PieRatioDataItem() { Title = "s2", Value = 30 },
                //new PieDataItem() { Title = "s3", Value = 20 },
                //new PieDataItem() { Title = "s4", Value = 80 }
            };

        int longestName = 0;
        private string[] titles;
        public void UpdatePeptideChart(double[][] intensities, string[][] scanTimes, string[] names)
        {
            Dictionary<string, double[]> intByScan = new Dictionary<string, double[]>();
            for (int i = 0; i < names.Length; i++)
            {
                for (int j = 0; j < scanTimes[i].Length; j++)
                {
                    if (!intByScan.ContainsKey(scanTimes[i][j]))
                        intByScan.Add(scanTimes[i][j], new double[names.Length]);
                    intByScan[scanTimes[i][j]][i] = intensities[i][j];
                }
                if (names[i].Length > longestName)
                    longestName = names[i].Length;
            }

            this.titles = names;

            foreach(string key in intByScan.Keys)
                data.Add(new PeptideDataItem() { Intensities = intByScan[key], ScanTime = key });
            _dataLoaded = false;
        }
                        //for (int j = 0; j < intensities[i].Length; j++)
                        //data.Add(new PeptideDataItem() { Title = names[i], Intensity = intensities[i][j], ScanTime = scanTimes[i][j] });                
            
         //       if (names[i].Length > longestName)
         //           longestName = names[i].Length;
         //   }
         //   this.titles = names;
         //   _dataLoaded = false;
        //}

        public SynthPeptideUC(double[][] intensities, string[][] scanTimes, string[] names, string categoryTitle)
        {
            UpdatePeptideChart(intensities, scanTimes, names);
            InitializeComponent();

            Title.Text = categoryTitle;
        }

        private bool _dataLoaded = false;
        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            if (!_dataLoaded)
            {
                serial1.Graphs.Clear();
                for(int i = 0; i < titles.Length; i++)
                {
                    AmCharts.Windows.QuickCharts.LineGraph line = new AmCharts.Windows.QuickCharts.LineGraph() { ValueMemberPath = "Intensities[" + i + "]", Title = titles[i] };
                    serial1.Graphs.Add(line);
                }
                serial1.DataSource = data;
                serial1.MinWidth = 600 + longestName * 6;
                _dataLoaded = true;
            }
        }
        /*
        public RatioUC()
        {
            InitializeComponent();            
            Title.Text = "Precursor";
        }//*/
    }

    public class PeptideDataItem : INotifyPropertyChanged
    {
        //public string Title { get; set; }
        //public string Category { get; set; }

        private double[] _Intensities;
        public double[] Intensities
        {
            get
            {
                return _Intensities;
            }
            set
            {
                _Intensities = value;

                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("Intensities"));
                }
            }
        }

        private string _ScanTime;
        public string ScanTime
        {
            get
            {
                return _ScanTime;
            }
            set
            {
                _ScanTime = value;

                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("ScanTime"));
                }
            }
        }

        #region INotifyPropertyChanged Members

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion
    }
}
