#ifndef KMP_HPP
#define KMP_HPP

template<typename DataRandomIterator, typename IndexRandomIterator, typename IndexT = typename std::iterator_traits<IndexRandomIterator>::value_type>
void BuildPartialMatchTable(DataRandomIterator data, IndexT count, IndexRandomIterator index){
  *index = (IndexT)-1;
  IndexT j = 0;
  IndexT k = (IndexT)-1;
  while(j<count-1){
    if(k == (IndexT)-1 || *(data+j)==*(data+k)){
      ++j;
      *(index+j) = ++k;
    }else
      k = *(index+k);
  }
}

template<typename DataRandomIterator, typename IndexT>
IndexT KMPSearch(
    DataRandomIterator data, IndexT data_count,
    DataRandomIterator pattern, IndexT pattern_count
                             ){
  IndexT* lps = new IndexT[pattern_count];
  BuildPartialMatchTable(pattern,pattern_count,lps);

  IndexT i=0, j=0;
  while(i<data_count){
    if(j==(IndexT)-1 || *(pattern+j)==*(data+i)){
      ++i;
      ++j;
      if(j==pattern_count)
        return i-j;
    }
    else{
      j = lps[j];
    }
  }
  return (IndexT)-1;
}

// const void* find(
//     const void* a, const void* b,
//     size_t size, size_t count_a, size_t count_b
//     , int (*compar)(const void*, const void*, size_t ) ){
//   const char* pa = (const char*)a;
//   const char* pb = (const char*)b;

//   size_t partial_match_index[count_b + 1];
//   memset(partial_match_index, 0, sizeof(size_t)*(count_b+1));

//   for (size_t i = 1; i < count_b; ++i)
//   {
//     size_t j = partial_match_index[i + 1];
//     while (j > 0 && compar(pb+size*j, pb+size*i, size)!=0) {
//       j = partial_match_index[j];
//     }
//     if (j > 0 || compar(pb+size*j, pb+size*i, size)==0) {
//       partial_match_index[i+1] = j+1;
//     }
//   }

//   size_t i = 0, j = 0;
//   while( i < count_a )
//   {
//     if(compar(pa+size*i, pb+size*j, size)==0)
//     {
//       if (++j == count_b) {
//         return pa + (i-j+1)*size;
//       }
//     }
//     else if (j > 0) {
//       j = partial_match_index[j];
//       continue;
//     }
//     ++i;
//   }
//   return nullptr;
// }

#endif
