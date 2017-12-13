from ROOT import TKey, TObject 

def FileIterator(d, basepath="/"):
    for key in d.GetListOfKeys():
        kname = key.GetName()

        if (key.IsFolder() and key.ReadObj().ClassName() != 'TTree'): 
            for i in FileIterator(d.Get(kname), basepath+kname+"/"):
                yield i
        else:
            yield basepath+kname, d.Get(kname)
        
                
