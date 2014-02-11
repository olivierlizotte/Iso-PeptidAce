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
    public partial class RatioUC : UserControl
    {
        ObservableCollection<PieRatioDataItem> data = new ObservableCollection<PieRatioDataItem>()
            {
                //new PieRatioDataItem() { Title = "s1", Value = 10 },
                //new PieRatioDataItem() { Title = "s2", Value = 30 },
                //new PieDataItem() { Title = "s3", Value = 20 },
                //new PieDataItem() { Title = "s4", Value = 80 }
            };

        int longestName = 0;
        public void UpdatePieChart(double[] ratios, string[] names)
        {
            for (int i = 0; i < ratios.Length; i++)
            {
                data.Add(new PieRatioDataItem() { Title = names[i], Value = ratios[i] });
                if (names[i].Length > longestName)
                    longestName = names[i].Length;
            }
        }

        public RatioUC(double[] ratios, string[] names, string name)
        {
            UpdatePieChart(ratios, names);
            InitializeComponent();

            Title.Text = name;
        }

        private void Window_Loaded(object sender, RoutedEventArgs e)
        {
            pie1.DataSource = data;            
            pie1.MinWidth = 200 +longestName * 6;
        }
        /*
        public RatioUC()
        {
            InitializeComponent();            
            Title.Text = "Precursor";
        }//*/
    }

    public class PieRatioDataItem : INotifyPropertyChanged
    {
        public string Title { get; set; }

        private double _value;
        public double Value
        {
            get
            {
                return _value;
            }
            set
            {
                _value = value;

                if (PropertyChanged != null)
                {
                    PropertyChanged(this, new PropertyChangedEventArgs("Value"));
                }
            }
        }

        #region INotifyPropertyChanged Members

        public event PropertyChangedEventHandler PropertyChanged;

        #endregion
    }
}
