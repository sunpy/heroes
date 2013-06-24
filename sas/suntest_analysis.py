import sas
import pandas
from datetime import datetime, timedelta
import matplotlib.pyplot as plt
import numpy as np

#filename = sas.create_summary_file(fits_files_dir)
#lc = sas.read_summary_file('eggs.csv')
lc = sas.read_summary_file('sas_summary_file.csv')

fig_filename = 'suntest_analysis_timeperiod'

timeperiod_num = 0

plot_num = 1
a = lc['CTLSOLUTION1'].plot()
a.set_title('CTL Solution X')
a.set_ylabel('Degrees')
plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
plt.show()

plot_num+=1
a = lc['CTLSOLUTION2'].plot()
a.set_title('CTL Solution X')
a.set_ylabel('Degrees')
plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
plt.show()

plot_num+=1
a = lc['SUN-CENTER1'].plot()
a.set_title('Sun Center X')
a.set_ylabel('Pixels')
plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
plt.show()

plot_num+=1
a = lc['SUN-CENTER2'].plot()
a.set_title('Sun Center Y')
a.set_ylabel('Pixels')
plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
plt.show()

for i in np.arange(1,4):
    if i == 1:
        start_time = datetime(2013,4,21,19,10)
        end_time = datetime(2013,4,21,19,30)
    if i == 2:
        start_time = datetime(2013,4,21,19,42)
        end_time = datetime(2013,4,21,20,12)
    if i == 3:
        start_time = datetime(2013,4,21,20,24)
        end_time = datetime(2013,4,21,20,45)

    plot_num = 0
    timeperiod_num = i
    
    plc=lc.ix[lc.index.indexer_between_time(start_time, end_time)]

    #convert from degrees to arcsec
    plc['CTLSOLUTION1'] = plc['CTLSOLUTION1']*60*60
    plc['CTLSOLUTION2'] = plc['CTLSOLUTION2']*60*60

    plot_num+=1
    l = plc['CTLSOLUTION1']
    ax = l.plot()
    #ax.set_ybound(-10,100)
    ax.set_ylabel('arcsec')
    ax.set_title('CTL Solution X period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    plot_num+=1
    ax = plc['CTLSOLUTION2'].plot()
    ax.set_ylabel('arcsec')
    #ax.set_ybound(2800,3200)
    ax.set_title('CTL Solution Y period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    plot_num+=1
    ax = plc['SUN-CENTER1'].plot()
    ax.set_ylabel('pixel')
    #ax.set_ybound(1020,1080)
    ax.set_title('Sun Center X period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    plot_num+=1
    ax = plc['SUN-CENTER2'].plot()
    ax.set_ylabel('pixel')
    #ax.set_ybound(425,440)
    ax.set_title('Sun Center Y period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    plot_num+=1
    r = np.sqrt(plc['CTLSOLUTION1']**2 + plc['CTLSOLUTION2']**2)
    ax = r.plot()
    ax.set_ylabel('arcsec')
    #ax.set_ybound(2800,3200)
    ax.set_title('CTL Solution period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    fit = np.polyfit(r.index.astype(np.int64), r.values,1)
    ylin = fit[0]*r.index.astype(np.int64) + fit[1]

    fit = pandas.TimeSeries(ylin, index=r.index)
    df = pandas.DataFrame(r)
    df[1] = fit
    df.columns = ['data', 'fit']

    plot_num+=1
    ax = df.plot()
    ax.set_ylabel('arcsec')
    #ax.set_ybound(2800,3200)
    ax.set_title('CTL Solution period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    diff = df['data'] - df['fit']

    plot_num+=1
    ax = diff.plot()
    ax.set_ylabel('arcsec')
    ax.set_ybound(-50,50)
    ax.set_title('Difference period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    plot_num+=1
    ax = diff.hist(bins=1000)
    ax.set_xlabel('arcsec')
    ax.set_xbound(-50,50)
    ax.set_title('Difference period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    plot_num+=1
    ax=plc['SUN-CENTER1'].plot()
    ax.set_ylabel('pixel')
    #ax.set_ybound(1000,1100)
    ax.set_title('Sun Center X period ' + str(timeperiod_num))
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()
    
    plot_num+=1
    ax=plc['SUN-CENTER2'].plot()
    ax.set_ylabel('pixel')
    ax.set_title('Sun Center Y period ' + str(timeperiod_num))
    #ax.set_ybound(420,450)
    plt.savefig(fig_filename + str(timeperiod_num) + '_' + str(plot_num))
    plt.show()

