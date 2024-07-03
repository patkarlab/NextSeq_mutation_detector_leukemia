import pandas as pd
import sys


args = sys.argv

filename = args[1]
outfile = args[2]


df = pd.read_csv(filename)
#cols=['Chr','Start','End','Ref','Alt','FILTER','SOMATIC FLAG','NUM_TOOLS','Tools.calledVariants','VAF','REF_COUNT','ALT_COUNT','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','Gene_full_name.refGene','Function_description.refGene','Disease_description.refGene','cosmic84','PopFreqMax','1000G_ALL','ExAC_ALL','CG46','ESP6500siv2_ALL','InterVar_automated']
common_cols=['Chr','Start','End','Ref','Alt','Variant_Callers','FILTER','SOMATIC_FLAG','VariantCaller_Count','REF_COUNT','ALT_COUNT','VAF','Func.refGene','Gene.refGene','ExonicFunc.refGene','AAChange.refGene','Gene_full_name.refGene','Function_description.refGene','Disease_description.refGene','cosmic84','PopFreqMax','1000G_ALL','ExAC_ALL','CG46','ESP6500siv2_ALL','InterVar_automated']
x = df['Otherinfo1']
discarded_column=df.columns.get_loc('Otherinfo2')
data = dict()
print('Filename:',filename,'Orginal shape',df.shape)
data.setdefault('FILTER', '-1')
data.setdefault('SOMATIC_FLAG', '-1')
data.setdefault('VariantCaller_Count',[])
data.setdefault('Variant_Callers',[])
data.setdefault('REF_COUNT', [])
data.setdefault('ALT_COUNT', [])
data.setdefault('VAF', [])
vartools =set()
for row in x:
    rowitems=row.split('\t')
    info=rowitems[10].split(';')
    set_str = [val for key,val in  enumerate(info) if val.startswith('set=')][0].split('=')[1]
    vartools.add(set_str)
    c = 0
    #check set values according put for other callers
    if(set_str =='variant' or set_str =='variant-filterInvariant2'):
        data['Variant_Callers'].append('Freebayes')
        #data['VariantCaller_Count'].append('1')
        c = 1
    elif (set_str =='Intersection'):
        data['Variant_Callers'].append('Freebayes | Platypus')
        #data['VariantCaller_Count'].append('2')
        c = 2
    else:
        data['Variant_Callers'].append('Platypus')
        #data['VariantCaller_Count'].append('1')
        c = 1
    formatval=rowitems[-1].split(':')
    refreads=formatval[-1]
    refcountlist=list(map(int,refreads.split(',')))
    refcount=sum(refcountlist)
    data['REF_COUNT'].append(refcount)
    formatindex=rowitems[-2].split(':')
    altreads=formatval[formatindex.index('AO')]
    altcountlist=list(map(int,altreads.split(',')))
    altcount=sum(altcountlist)
    data['ALT_COUNT'].append(altcount)
    allele_fraction=float(int(altcount)/int(altcount+refcount))
    vaf="{:.2%}".format(allele_fraction)
    data['VAF'].append(vaf)
    data['VariantCaller_Count'].append(int(c))
print('set values:',vartools)
df1=df.iloc[:,:5]
df2=pd.DataFrame(data, columns=data.keys())
df3=df.iloc[:,5:discarded_column]
horizontal_stack = pd.concat([df1, df2, df3], axis=1)
# Format changes
horizontal_stack['cosmic84']=horizontal_stack['cosmic84'].astype(str).str.replace(',' , ';')
horizontal_stack['AAChange.refGene']=horizontal_stack['AAChange.refGene'].str.replace(',' , ';')
horizontal_stack.replace(to_replace='.', value='-1', inplace=True)
horizontal_stack=horizontal_stack.reindex(columns = common_cols)
#horizontal_stack.rename(columns = {'Func.refGene':'Variant_Site', 'ExonicFunc.refGene':'Variant_Function', 'AAChange.refGene':'Amino.Acid_Change','Gene.refGene':'Gene','Gene_full_name.refGene':'Gene.Full_Name','Function_description.refGene':'Gene.Function','Disease_description.refGene':'Disease_description'}, inplace = True)
#Data Filtering..
#df=df.loc[(df.Variant_Site == 'exonic') & (df.Variant_Function != 'synonymous SNV')]
print('New shape:',horizontal_stack.shape)
horizontal_stack.to_csv(outfile, index=False)
# TODO Comment this
#horizontal_stack.to_csv(filepath+'OCIAMl3.combined.csv', index=False)
