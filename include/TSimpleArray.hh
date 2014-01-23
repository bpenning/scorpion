#ifndef TSimpleArray_h
#define TSimpleArray_h

#include "TObjArray.h"
#include "TRefArray.h"

#include <map>

class TCompare;

class TBaseArray
{
  public:
  TBaseArray(Int_t size = 100) :
    fIsOwner(kTRUE), fIsRefArray(kFALSE), fArray(0), fIter(0)
  {
    fArray = new TObjArray(size);
    fIter = fArray->MakeIterator();
  }
  
  TBaseArray(TObjArray *array) :
    fIsOwner(kFALSE), fIsRefArray(kFALSE), fArray(array), fIter(0)
  {
    fIter = fArray->MakeIterator();
  }

  TBaseArray(TRefArray &array) :
    fIsOwner(kFALSE), fIsRefArray(kTRUE), fArray(&array), fIter(0)
  {
    fIter = fArray->MakeIterator();
  }

  virtual ~TBaseArray()
  {
    delete fIter;
    if(fIsOwner) delete fArray;
  }

  Int_t GetEntries() const
  {
    if(fIsRefArray)
    {
      return static_cast<TRefArray*>(fArray)->GetEntriesFast();
    }
    else
    {
      return static_cast<TObjArray*>(fArray)->GetEntriesFast();
    }
  }    

  void Clear()
  {
    fArray->Clear();
    Reset();
    ClearFilter();
  }

  void Reset() { fIter->Reset(); }

  virtual void ClearFilter(const void *filter = 0) = 0;

private:

  Bool_t fIsOwner;

protected:

  Bool_t fIsRefArray;
  TSeqCollection *fArray;
  TIterator *fIter;
};

//------------------------------------------------------------------------------

template<typename T>
class TSimpleArray: public TBaseArray
{
public:

  class FilterExeption{};

  TSimpleArray(Int_t size = 100) :
    TBaseArray(size), fMap() {}

  TSimpleArray(TObjArray *array) :
    TBaseArray(array), fMap() {}

  TSimpleArray(TRefArray &array) :
    TBaseArray(array), fMap() {}
  
  ~TSimpleArray()
  {
    TFilterMapConstIter it_map;
    TCategoryMapConstIter it_submap;
    for(it_map = fMap.begin(); it_map != fMap.end(); ++it_map)
    {
      for(it_submap = it_map->second.second.begin();
          it_submap != it_map->second.second.end(); ++it_submap)
      {
        delete (it_submap->second);
      }
    }
  }

  void ClearFilter(const void *filter = 0)
  {
    TFilterMapIter it_map;
    TCategoryMapIter it_submap;
    if(filter)
    {
      it_map = fMap.find(filter);
      if(it_map != fMap.end())
      {
        it_map->second.first = kTRUE;
        for(it_submap = it_map->second.second.begin();
            it_submap != it_map->second.second.end(); ++it_submap)
        {
          it_submap->second->Clear();
        }
      }
    }
    else
    {
      for(it_map = fMap.begin(); it_map != fMap.end(); ++it_map)
      {
        it_map->second.first = kTRUE;
        for(it_submap = it_map->second.second.begin();
            it_submap != it_map->second.second.end(); ++it_submap)
        {
          it_submap->second->Clear();
        }
      }
    }
  }
  
  void Add(T *obj) { fArray->Add(obj); }

  T *At(Int_t i) const
  {
    if(fIsRefArray)
    {
      return static_cast<T*>((*static_cast<TRefArray*>(fArray))[i]); 
    }
    else
    {
      return static_cast<T*>((*static_cast<TObjArray*>(fArray))[i]); 
    }
  }
  T *operator[](Int_t i) const { return At(i); }

  void Sort(TCompare *comp = 0)
  {
    TCompare *temp = T::fgCompare;
    if(comp)
    {
      T::fgCompare = comp;
    }
    if(fIsRefArray)
    {
      static_cast<TRefArray*>(fArray)->Sort();
    }
    else
    {
      static_cast<TObjArray*>(fArray)->Sort();
    }
    if(comp)
    {
      T::fgCompare = temp;
    }
  }

  template <typename F>
  TSimpleArray<T> *GetSubArray(const F *filter, Int_t category = 1)
  {
    T *element;
    Int_t result;
    TCategoryMapIter it_submap;
    std::pair<TCategoryMapIter, bool> pair_submap;
    std::pair<TFilterMapIter, bool> pair_map;

    TFilterMapIter it_map = fMap.find(filter);
    if(it_map == fMap.end())
    {
      pair_map = fMap.insert(make_pair(filter, make_pair(kTRUE, TCategoryMap())));
      if(!pair_map.second) throw FilterExeption();

      it_map = pair_map.first;
    }

    if(it_map->second.first)
    {
      it_map->second.first = kFALSE;
      Reset();
      while(element = Next())
      {
        result = (*filter)(element);
        if(result < 0) continue;
        it_submap = it_map->second.second.find(result);
        if(it_submap == it_map->second.second.end())
        {
          pair_submap = it_map->second.second.insert(make_pair(result,
                                                               new TSimpleArray<T>(fArray->GetSize())));
          if(!pair_submap.second) throw FilterExeption();

          pair_submap.first->second->Add(element);
        }
        else
        {
          it_submap->second->Add(element);
        }
      }
    }

    it_submap = it_map->second.second.find(category);
    return (it_submap != it_map->second.second.end()) ? it_submap->second : 0;
  }

  T *Next() { return static_cast<T*>(fIter->Next()); }

private:
  typedef std::map<Int_t, TSimpleArray<T>*> TCategoryMap;
  typedef std::map< const void*, std::pair<Bool_t, TCategoryMap> > TFilterMap;

  typedef typename TFilterMap::iterator TFilterMapIter;
  typedef typename TCategoryMap::iterator TCategoryMapIter;
  typedef typename TFilterMap::const_iterator TFilterMapConstIter;
  typedef typename TCategoryMap::const_iterator TCategoryMapConstIter;

  TFilterMap fMap;
};

#endif /* TSimpleArray */

