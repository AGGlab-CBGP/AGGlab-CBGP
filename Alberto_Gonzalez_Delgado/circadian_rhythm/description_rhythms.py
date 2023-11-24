#!/bin/python

import sys
import argparse
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin
from sklearn.pipeline import make_pipeline
from sklearn.linear_model import LinearRegression
from scipy.signal import find_peaks
from scipy.stats import skew, kurtosis 

def check_arg(args=None):
    '''
	Description:
        Function collect arguments from command line using argparse
    Input:
        args # command line arguments
    Constant:
        None
    Variables
        parser
    Return
        parser.parse_args() # Parsed arguments
    '''
    parser = argparse.ArgumentParser(prog = 'oscillatory_analysis.py', formatter_class=argparse.RawDescriptionHelpFormatter, description= 'oscillatory_analysis.py analyses the power spectra of the raw data for detecting oscillatory pattern')
    parser.add_argument('--sample','-s',required=True,help='Insert name of the sample (lightconditions_genotype)')
    return parser.parse_args()

class GaussianFeatures(BaseEstimator, TransformerMixin):
    """Uniformly spaced Gaussian features for one-dimensional input"""
    def __init__(self, N, width_factor):
        self.N = N
        self.width_factor = width_factor
    @staticmethod
    def _gauss_basis(x, y, width, axis=None):
        arg = (x - y) / width
        return np.exp(-0.5 * np.sum(arg ** 2, axis))
    def fit(self, X, y=None):
        # create N centers spread along the data range
        self.centers_ = np.linspace(X.min(), X.max(), self.N) #Aqui esta el error, no está accediendo al xmin y xmax
        if len(self.centers_) < 2:
            raise ValueError("Number of centers is too small")
        self.width_ = self.width_factor * (self.centers_[1] - self.centers_[0])
        return self
    def transform(self, X):
        return self._gauss_basis(X[:, :, np.newaxis], self.centers_,
                                 self.width_, axis=1)
## Data_set
class data_set:
    def __init__(self,condition,sep):
        self.condition=condition
        self.sep=sep
    def read_raw_file(self):
        '''
            Description:
                Function normalized counts raw expression file
            Input:
                self
            Return:
                data frame
        '''
        self.raw=pd.read_csv(f'./exp/{self.condition}.csv',sep=self.sep)
        return self.raw.iloc[:, 1:]
    def read_time_file(self):
        '''
            Description:
                Function read time  file
            Input:
                self
            Return:
                data frame
        '''
        self.time=pd.read_csv(f'./time/{self.condition}_time.csv',sep=self.sep)
        return self.time.iloc[:, 1:]
    def proc_file(self):
        '''
            Description:
                Function processed data file
            Input:
                self
            Return:
                data frame
        '''
        self.proc=pd.read_csv(f'proc_data/{self.condition}.csv',sep=self.sep)
        return self.proc.iloc[:, 1:]

    def read_normalized_file(self):
        '''
            Description:
                Function processed data file
            Input:
                self
            Return:
                data frame
        '''
        self.proc=pd.read_csv(f'normalized_expression/{self.condition}_normalized.csv',sep=self.sep)
        return self.proc.iloc[:, 1:]
        
    def preprocess(self,raw,time_file,time_interp,N,width_factor,write=False):
            '''
                Description:
                    Function normalizes gene expression. It read the parsed raw expression file
                Input:
                    self , time (int), time_interp (float),gene_list (array), write (boolean)
                Return:
                    data frame and csv file (optional)
            '''
            residuals = pd.DataFrame(columns=['gene', 'centers', 'coef'])
            dict_={}
            xfit=time_interp #np.arange(0,24,0.166667)
            gene_list=raw.columns.values
            for gene in gene_list:
                rep1_zipped=zip(raw[gene].values[:12],time_file[gene].values)
                rep1_zipped_sorted = sorted(rep1_zipped, key=lambda x: x[1])
                rep1=[x[0] for x in rep1_zipped_sorted]
                rep2_zipped=zip(raw[gene].values[12:],time_file[gene].values)
                rep2_zipped_sorted = sorted(rep2_zipped, key=lambda x: x[1])
                rep2=[x[0] for x in rep2_zipped_sorted]
                x=[x[1] for x in rep2_zipped_sorted]
                x_array = np.array(x)[:, np.newaxis]
                zipped=zip(rep1,rep2)
                means = [np.mean(x) for x in zipped]
                gauss_model = make_pipeline(GaussianFeatures(N,width_factor),LinearRegression())
                gauss_model.fit(x_array, means)
                exp_new = gauss_model.predict(xfit[:, np.newaxis])
                dict_[gene]=exp_new
                self.preprocessed=pd.DataFrame(dict_)
                gene=gene
                centers=gauss_model.steps[0][1].centers_
                coef=gauss_model.steps[1][1].coef_
                temp=[gene,pd.Series(centers),pd.Series(coef)]
                residuals.loc[len(residuals)] = temp
            if write:
                pd.DataFrame(dict_).to_csv(f'proc_data/{self.condition}.csv')
                residuals.to_csv(f'proc_data/{self.condition}_residuals.csv')
            return self.preprocessed , residuals

    def normalize(self, raw_exp, time_file, read=True, write=False):
        '''
        Description:
            Function normalizes gene expression. It can read the processed expression file or call the processing function
            to process data
        Input:
            self , gene_list (array), read (boolean), write (boolean)
        Return:
            data frame and csv file (optional) 
        ''' 
        if read:
            prep = self.proc_file()
        else:
            prep = self.preprocess(raw_exp, time_file,np.arange(1, 24, 0.1666666667))
        dict = {}
        gene_list=prep.columns.values
        for gene in gene_list:
            expression = prep[gene]
            min_value = np.min(expression)
            max_value = np.max(expression)
            normalized_exp = (expression - min_value) / (max_value - min_value)
            dict[gene] = normalized_exp.tolist()

        if write:
            pd.DataFrame(dict).to_csv(f'normalized_expression/{self.condition}_normalized.csv', index=False)

        self.normalized = pd.DataFrame(dict)
        return self.normalized
# Oscillation description
class oscillatory_description:
    
    def __init__(self,condition):
        self.condition=condition

    def read_peliminary_peaks_descriptors(self):
        '''
    	    Description:
                Function read preliminary descriptors file
            Input:
                self
            Return:
                data frame 
        ''' 
        self.preliminary_descriptors=pd.read_csv(f'preliminary_peaks/peaks_df_{self.condition}.csv',sep=',')
        return self.preliminary_descriptors.iloc[:, 1:]
    
    def read_descriptors(self):
        '''
    	    Description:
                Function read descriptors file
            Input:
                self
            Return:
                data frame 
        ''' 
        self.descriptors=pd.read_csv(f'peaks/features_{self.condition}.csv')
        return self.descriptors

    def finding_peaks(self):
        '''
    	    Description:
                Function detects peaks from gene expression
            Input:
                self, gene list (string) 
            Return:
                dictionary
        ''' 
        peaks_array={}
        expression=data_set(f'{self.condition}',',').proc_file()
        gene_list=expression.columns.values
        for gene in gene_list:
            peaks_array[gene]=find_peaks(expression[gene],prominence=0,threshold=0,width=0,height=0,plateau_size=0)
        self.peak_a=peaks_array
        return self.peak_a
    
    def calculate_preliminary_descriptors(self,write=False):
        '''
    	    Description:
                Function calculates preliminary raw descriptors from detected peaks
            Input:
                self, gene list (string) and write (boolean) to indicate wether the output file is generated
            Return:
                data frame or output.csv file (optional)
        '''  
        peaks=[];prominences=[];heights=[];plateau_size=[];widths=[];peaks_count=[];genes=[];right_t=[];left_t=[]
        peak_a=self.finding_peaks()
        expression=data_set(f'{self.condition}',',').proc_file()
        gene_list=expression.columns.values
        for i in range(len(gene_list)):
            genes.append(gene_list[i])
            peaks.append(peak_a[gene_list[i]][0])
            peaks_count.append(len(peak_a[gene_list[i]][0]))
            prominences.append(peak_a[gene_list[i]][1]['prominences'])
            heights.append(peak_a[gene_list[i]][1]['peak_heights'])
            plateau_size.append(peak_a[gene_list[i]][1]['plateau_sizes'])
            widths.append(peak_a[gene_list[i]][1]['widths'])
            right_t.append(peak_a[gene_list[i]][1]['right_ips'])
            left_t.append(peak_a[gene_list[i]][1]['left_ips'])
        features={
            'genes' : genes,
            'peaks' : peaks,
            'peaks_count' : peaks_count,
            'prominences' : prominences,
            'heights' : heights,
            'plateau_size' : plateau_size,
            'widths' : widths,
            'time_peak_start' : left_t,
            'time_peak_finishes' : right_t
        }
        self.preliminary_descriptors=pd.DataFrame(features)
        if write:
            self.preliminary_descriptors.to_csv(f'preliminary_peaks/peaks_df_{self.condition}.csv')
        return self.preliminary_descriptors

    
    def calculate_descriptors(self,write=False):
        '''
    	    Description:
                Function extract features and convert them from string into numeric and calculate other features
            Input:
                self, gene list (string) and write (boolean) to indicate wether the output file is generated
            Return:
                data frame or output.csv file (optional)
        '''  
        self.preliminary_descriptors=self.read_peliminary_peaks_descriptors()
        self.exp_normalized=data_set(self.condition,',').read_normalized_file()
        gene_list=self.exp_normalized.columns.values
        self.exp_prep=data_set(self.condition,',').proc_file()
        residuals=pd.read_csv(f'./proc_data/{self.condition}_residuals.csv',sep=',')
        phases_a=[];prominences_a=[];plateau_size_a=[];widths_a=[];peaks_count_a=[];gene_a=[];increase_a=[];start_a=[];end_a=[];heights_a=[];auc_a=[];p_value_a=[];phases_a=[];peak_increasing=[];peak_decreasing=[];phase_oscillation_a=[];amplitude_a=[];period_a=[];auc_exp_a=[]
        centers_a=[];coefs_a=[];skew_a=[];kurt_a=[]
        time=np.arange(1,24,0.16666667)
        for i in range(len(self.preliminary_descriptors.genes)):
            gene=self.preliminary_descriptors.genes[i]
            if self.preliminary_descriptors.genes[i] in gene_list:
              #Appending features from the other table generatged previously
    
              ######################################### 0 peaks #################################################
    
              #Eliminating housekeeping genes from the analysis
              if self.preliminary_descriptors.peaks_count[i]==0:
                  continue
              
              ######################################### >1 peaks #################################################
              
              elif self.preliminary_descriptors.peaks_count[i]>0:
                  #Cleaning input data and extracting peaks and widths. Calculating mean of phase(peaks),widhts,prominences,heights,plateau_size.
                  peaks=self.preliminary_descriptors.peaks[i].split(" ")
                  widths=self.preliminary_descriptors.widths[i].split(" ")
                  prominences=self.preliminary_descriptors.prominences[i].split(" ")
                  plateau_size=self.preliminary_descriptors.plateau_size[i].split(" ")
                  height=self.preliminary_descriptors.heights[i].split(" ")
                  start=self.preliminary_descriptors.time_peak_start[i].split(" ")
                  final=self.preliminary_descriptors.time_peak_finishes[i].split(" ")
                  expression=self.exp_normalized[self.preliminary_descriptors.genes[i]]
                  exp_no_normalized=self.exp_prep[self.preliminary_descriptors.genes[i]]
                  centers=residuals[residuals['gene']==gene]['centers'].values[0].replace('dtype: float64','').replace('\n','').split()[1:]
                  centers=list(map(float,centers))
                  coefs=residuals[residuals['gene']==gene]['coef'].values[0].replace('dtype: float64','').replace('\n','').split()[1:]
                  coefs=list(map(float,coefs))
                  peaks_c=[]
                  for item in peaks:
                      try:
                          peaks_c.append(float(item.replace("[","").replace("]","").replace(" ","")))
                      except:
                          continue
                  widths_c=[]
                  for item in widths:
                      try:
                          widths_c.append(float(item.replace("[","").replace("]","").replace(" ","")))
                      except:
                          continue
                      
                  prominences_c=[]
                  for item in prominences:
                      try:
                          prominences_c.append(float(item.replace("[","").replace("]","").replace(" ","")))
                      except:
                          continue
                      
                  plateau_size_c=[]
                  for item in plateau_size:
                      try:
                          plateau_size_c.append(float(item.replace("[","").replace("]","").replace(" ","")))
                      except:
                          continue
                      
                  height_c=[]
                  for item in height:
                      try:
                          height_c.append(float(item.replace("[","").replace("]","").replace(" ","")))
                      except:
                          continue
                      
                  start_c=[]
                  for item in start:
                      try:
                          start_c.append(float(item.replace("[","").replace("]","").replace(" ","")))
                      except:
                          continue
                      
                  final_c=[]
                  for item in final:
                      try:
                          final_c.append(float(item.replace("[","").replace("]","").replace(" ","")))
                      except:
                          continue
                      
                  peaks_processed=[]
                  widths_processed=[]
                  prominences_processed=[]
                  plateau_size_processed=[]
                  height_processed=[]
                  start_processed=[]
                  final_processed=[]
                  symmetry=[]
                  
                  highest_prominence=max(prominences_c)
                  for p in range(0,len(peaks_c)):
                     if prominences_c[p]==highest_prominence:
                        if (peaks_c[p]*0.1666666667+1)<24:
                          peaks_processed.append((peaks_c[p]*0.1666666667)+1)
                        else:
                           peaks_processed.append((peaks_c[p]*0.1666666667)-24+1) 
                        start_processed.append(start_c[p]*0.1666666667+1)
                        final_processed.append(final_c[p]*0.1666666667+1)
                        widths_processed.append(widths_c[p]*0.1666666667)
                        prominences_processed.append(prominences_c[p])
                        plateau_size_processed.append(plateau_size_c[p])
                        height_processed.append(height_c[p])
                        right=round(np.mean(final_c[p]*0.1666666667)-np.mean(peaks_c[p]*0.1666666667),3)
                        left=round(np.mean(peaks_c[p]*0.1666666667)-np.mean(start_c[p]*0.1666666667),3)
                        if right>left:
                        #Faster increase = 0
                        #Faster decrease = 1                   
                          symmetry.append(0)
                        else:
                          symmetry.append(1)
                        AUC_shape=np.trapz(expression[int(start_c[p]):int(final_c[p])])
                        AUC_exp=np.trapz(exp_no_normalized[int(exp_no_normalized[p]):int(final_c[p])])
                        distribution=[]
                        for i,t in enumerate(time):
                            distribution.extend([time[i]]*round(exp_no_normalized[i]*100))
                        skewness=skew(distribution, bias=True)
                        krt=kurtosis(distribution,bias=True)
                  phases=np.mean(peaks_processed)
                  phases_a.append(round(phases))
                  peak_increasing.append(left)
                  peak_decreasing.append(right)
                  witdths_m=np.mean(widths_processed)
                  widths_a.append(round(witdths_m,4))
                  #prominences_m=np.mean(prominences_processed)
                  #prominences_a.append(round(prominences_m,4))
                  prominences_a.append(np.mean(prominences_processed))
    
                  plateau_size_m=np.mean(plateau_size_processed)
                  plateau_size_a.append(round(plateau_size_m,4))
    
                  #heights_m=np.mean(height_processed)    
                  #heights_a.append(round(heights_m,2))
                 
                  gene_a.append(self.preliminary_descriptors.genes[i])
                  peaks_count_a.append(len(peaks_c))
                  heights_a.append(np.mean(height_processed))
                  increase_a.append(round(np.mean(symmetry),3))
                  start_a.append(np.mean(start_processed))
                  end_a.append(np.mean(final_processed))
                  auc_a.append(AUC_shape)
                  auc_exp_a.append(AUC_exp)
                  skew_a.append(skewness)
                  kurt_a.append(krt)
                  centers_a.append(centers)
                  coefs_a.append(coefs)
            else:
               continue
        features={
            'gene' : gene_a,
            'peaks_count' : peaks_count_a,
            'prominences' : prominences_a,
            'widths' : widths_a,
            'heights' : heights_a,
            'plateau_size' : plateau_size_a,
            'AUC_shape' : auc_a,
            'AUC_exp' : auc_exp_a,
            'time_increasing' : peak_increasing,
            'time_decreasing' : peak_decreasing,
            'centers' : centers_a,
            'coefs' : coefs_a,
            'skewness' : skew_a,
            'kurtosis' : kurt_a,
            'start' : start_a,
            'end' : end_a

        }
        if write:
            pd.DataFrame(features).to_csv(f'peaks/features_{self.condition}.csv')
        self.descriptors=pd.DataFrame(features)
        return self.descriptors
            

################################################################  Variables  #########################################################################################

arguments = check_arg(sys.argv[1:])
sample=arguments.sample

#################################################################  Analysis  #########################################################################################

raw_exp = data_set(sample, ',').read_raw_file()
time_file=data_set(sample,',').read_time_file()
prep_exp = data_set(sample, ',').preprocess(raw_exp, time_file,np.arange(1, 24, 0.1666666667),N=4,width_factor=10, write=True)
prep=data_set(sample,',').proc_file()
normalized=data_set(sample,',').normalize(raw_exp,time_file,read=True,write=True)
normalized=data_set(sample,',').read_normalized_file()
preliminary_descriptors=oscillatory_description(sample).calculate_preliminary_descriptors(write=True)
descriptors=oscillatory_description(sample).calculate_descriptors(write=True)

