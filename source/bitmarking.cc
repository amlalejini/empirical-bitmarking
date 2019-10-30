#include <iostream>
#include <chrono>
#include <array>
#include <tuple>
#include <map>
#include <list>

#include "base/vector.h"

#include "tools/BitSet.h"
#include "tools/BitVector.h"

#include "data/DataFile.h"

#include "tools/Random.h"
#include "tools/random_utils.h"

#include <thread>

// Replicates - 100

// Tests - Sizes
// - 16 Bit
// - 32 Bit
// - 64 Bit
// - 128 Bit
// - 256 Bit

// Tests - Logic operations
// - ==
// - <
// - >
// - GetUInt
// - CountOnes_Mixed
// - NOT
// - AND
// - OR
// - NAND
// - AND
// - NOR
// - XOR
// - EQU
// - SHIFT

size_t SEED = 2;
size_t NUM_REPLICATES = 100;    // How many replicates should we do?
size_t NUM_ITERATIONS = 10000;  // How many times to repeat each operation before measuring time elapsed?


class Benchmark {
public:


private:
  size_t num_iterations;
  size_t num_replicates;

  struct Info {
    size_t replicate=0;
    std::string treatment="";
    size_t bits=0;
    std::string operation="";
    double time=0.0;
    size_t n=0;
    double value=0.0;
    std::string container="";
  } info;

  emp::DataFile datafile;
  emp::Random random;

  emp::vector< std::tuple<emp::BitVector,emp::BitVector> > bit_vectors;
  emp::vector< std::tuple<emp::BitSet<16>,emp::BitSet<16>> > bit_set_16;
  emp::vector< std::tuple<emp::BitSet<31>,emp::BitSet<31>> > bit_set_31;
  emp::vector< std::tuple<emp::BitSet<32>,emp::BitSet<32>> > bit_set_32;
  emp::vector< std::tuple<emp::BitSet<64>,emp::BitSet<64>> > bit_set_64;
  emp::vector< std::tuple<emp::BitSet<128>,emp::BitSet<128>> > bit_set_128;
  emp::vector< std::tuple<emp::BitSet<256>,emp::BitSet<256>> > bit_set_256;
  emp::vector< std::tuple<emp::BitSet<10000>,emp::BitSet<10000>> > bit_set_10000;

  std::list< std::tuple<emp::BitVector,emp::BitVector> > list_of_bit_vectors;
  std::list< std::tuple<emp::BitSet<64>,emp::BitSet<64>> > list_of_bit_set_64;

  std::tuple<emp::BitVector,emp::BitVector> lone_bit_vector;
  std::tuple<emp::BitSet<64>,emp::BitSet<64>> lone_bit_set_64;


  const std::string container = "T";
  const std::string container_vectorized = "emp::vector<T>";
  const std::string container_list = "std::list<T>";

  // Test bitvectors without vectorization
  template<size_t NUM_BITS>
  void BenchmarkBitVectors() {

    info.container = container;

    // Resize all bit vectors
    const size_t width = NUM_BITS;
    info.bits = width;

    std::get<0>(lone_bit_vector).Resize(width);
    std::get<1>(lone_bit_vector).Resize(width);


    // Compute timings for each replicate.
    for (size_t rep = 0; rep < num_replicates; ++rep) {
      info.replicate = rep;
      // Do bit vector:
      info.treatment = "bit_vector";
      // (1) Randomize bit_vectors

      emp::RandomizeBitVector(std::get<0>(lone_bit_vector), random);
      emp::RandomizeBitVector(std::get<1>(lone_bit_vector), random);


      std::cout << "bitvector: " << std::get<0>(lone_bit_vector) << ' ' << std::get<1>(lone_bit_vector) << std::endl;;

      emp::BitVector recipient;

      // (2) Time each operation!

      // ------------ CountOnes_Mixed ------------
      info.operation = "CountOnes_Mixed";

      double comulative = 0;

      auto start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        emp::BitVector & bv1 = std::get<1>(lone_bit_vector);


        bv[i%bv.size()] = !bv[i%bv.size()];
        comulative += bv.CountOnes_Mixed();
      }
      auto end_time = std::chrono::high_resolution_clock::now();

      std::cout << comulative << std::endl;

      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();


      // ------------ NOT ------------
      info.operation = "NOT";

      start_time = std::chrono::high_resolution_clock::now();
      
      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        bv[i%bv.size()] = !bv[i%bv.size()];
        bv = bv.NOT();
         }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
      
      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;


      // ------------ AND ------------
      info.operation = "AND";
      start_time = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv0 = std::get<0>(lone_bit_vector);
        emp::BitVector & bv1 = std::get<1>(lone_bit_vector);
        
        bv0[i%bv0.size()] = !bv0[i%bv0.size()];


        bv0 = bv0.AND(bv1);
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
      
      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;
//      std::cout << "AND: " << recipient << std::endl;;


      // ------------ OR ------------
      info.operation = "OR";
      start_time = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv0 = std::get<0>(lone_bit_vector);
        emp::BitVector & bv1 = std::get<1>(lone_bit_vector);
        
        bv0[i%bv0.size()] = !bv0[i%bv0.size()];


        bv0 = bv0.OR(bv1);
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
      
      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;
    }
  }
  // Test bit sets without vectorization
  template<size_t NUM_BITS>
  void BenchmarkBitSets(std::tuple<emp::BitSet<NUM_BITS>,emp::BitSet<NUM_BITS>> & bit_set_tuple) {

    info.container = container;

    // Resize all bit vectors
    const size_t width = NUM_BITS;
    info.bits = width;

    // Compute timings for each replicate.
    for (size_t rep = 0; rep < num_replicates; ++rep) {
      info.replicate = rep;
      // Do bit set:
      info.treatment = "bit_set";
      // (1) Randomize bit_sets

      std::get<0>(bit_set_tuple).Randomize(random);
      std::get<1>(bit_set_tuple).Randomize(random);

      emp::BitSet<NUM_BITS> recipient;
      
      // ------------ CountOnes_Mixed ------------
      info.operation = "CountOnes_Mixed";

      double comulative = 0;

      auto start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_set_tuple);

        bs.Toggle(i%bs.size());
        comulative += bs.CountOnes_Mixed();
      }
      auto end_time = std::chrono::high_resolution_clock::now();

      std::cout << comulative << std::endl;

      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ NOT ------------
      info.operation = "NOT";
      start_time = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(bit_set_tuple);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(bit_set_tuple);
        
        bs0.Toggle(i%bs0.size());

        bs0 = bs0.NOT();

      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
  
      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_set_tuple);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

      // ------------ AND ------------
      info.operation = "AND";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(bit_set_tuple);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(bit_set_tuple);
        
        bs0.Toggle(i%bs0.size());

        bs0 = bs0.AND(bs1);
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_set_tuple);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

      // ------------ OR ------------
      info.operation = "OR";
      start_time = std::chrono::high_resolution_clock::now();
      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(bit_set_tuple);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(bit_set_tuple);
        
        bs0.Toggle(i%bs0.size());

        bs0 = bs0.OR(bs1);
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_set_tuple);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;
      
    }
  }
  
  // Test bitvectors with vectorization
  template<size_t NUM_BITS>
  void BenchmarkBitVectors_vectorized() {

    info.container = container_vectorized;

    // Resize all bit vectors
    const size_t width = NUM_BITS;
    info.bits = width;
    for (size_t i = 0; i < NUM_ITERATIONS; ++i) {
      std::get<0>(bit_vectors[i]).Resize(width);
      std::get<1>(bit_vectors[i]).Resize(width);
    }

    // Compute timings for each replicate.
    for (size_t rep = 0; rep < num_replicates; ++rep) {
      info.replicate = rep;
      // Do bit vector:
      info.treatment = "bit_vector";
      // (1) Randomize bit_vectors
      for (size_t i = 0; i < num_iterations; ++i) {
        emp::RandomizeBitVector(std::get<0>(bit_vectors[i]), random);
        emp::RandomizeBitVector(std::get<1>(bit_vectors[i]), random);
      }
      emp::BitVector recipient;

      // (2) Time each operation!

      // ------------ CountOnes_Mixed ------------
      info.operation = "CountOnes_Mixed";

      double comulative = 0;

      auto start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(bit_vectors[i]);
        bv[i%bv.size()] = !bv[i%bv.size()];
        comulative += bv.CountOnes_Mixed();
      }
      auto end_time = std::chrono::high_resolution_clock::now();

      std::cout << comulative << std::endl;

      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ NOT ------------
      info.operation = "NOT";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_vectors.size(); ++i) {
        emp::BitVector & bv = std::get<0>(bit_vectors[i]);
        bv[i%bv.size()] = !bv[i%bv.size()];
        bv = bv.NOT();
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
      
      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

      // ------------ AND ------------
      info.operation = "AND";
      start_time = std::chrono::high_resolution_clock::now();
   
      for (size_t i = 0; i < bit_vectors.size(); ++i) {
        emp::BitVector & bv0 = std::get<0>(bit_vectors[i]);
        emp::BitVector & bv1 = std::get<0>(bit_vectors[i]);

        bv0[i%bv0.size()] = !bv0[i%bv0.size()];
        bv0 = bv0.AND(bv1);
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

      // ------------ OR ------------
      info.operation = "OR";
      start_time = std::chrono::high_resolution_clock::now();
      
      for (size_t i = 0; i < bit_vectors.size(); ++i) {
        emp::BitVector & bv0 = std::get<0>(bit_vectors[i]);
        emp::BitVector & bv1 = std::get<0>(bit_vectors[i]);

        bv0[i%bv0.size()] = !bv0[i%bv0.size()];
        bv0 = bv0.OR(bv1);
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitVector & bv = std::get<0>(lone_bit_vector);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

    }
  }
  // Test bit sets with vectorization
  template<size_t NUM_BITS>
  void BenchmarkBitSets_vectorized(emp::vector< std::tuple<emp::BitSet<NUM_BITS>,emp::BitSet<NUM_BITS>> > & bit_sets) {
    
    info.container = container_vectorized;
    
    // Resize all bit vectors
    const size_t width = NUM_BITS;
    info.bits = width;

    // Compute timings for each replicate.
    for (size_t rep = 0; rep < num_replicates; ++rep) {
      info.replicate = rep;
      // Do bit set:
      info.treatment = "bit_set";
      // (1) Randomize bit_sets
      for (size_t i = 0; i < num_iterations; ++i) {
        std::get<0>(bit_sets[i]).Randomize(random);
        std::get<1>(bit_sets[i]).Randomize(random);
      }
      emp::BitSet<NUM_BITS> recipient;

      info.operation = "CountOnes_Mixed";

      double comulative = 0;

      auto start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_sets[i]);

        bs.Toggle(i%bs.size());
        comulative += bs.CountOnes_Mixed();
      }
      auto end_time = std::chrono::high_resolution_clock::now();

      std::cout << comulative << std::endl;

      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ NOT ------------
      info.operation = "NOT";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_sets.size(); ++i) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(bit_sets[i]);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(bit_sets[i]);
        
        bs0.Toggle(i%bs0.size());

        bs0 = bs0.NOT();
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_sets[i]);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

      // ------------ AND ------------
      info.operation = "AND";
      start_time = std::chrono::high_resolution_clock::now();

       for (size_t i = 0; i < bit_sets.size(); ++i) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(bit_sets[i]);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(bit_sets[i]);
        
        bs0.Toggle(i%bs0.size());

        bs0 = bs0.AND(bs1);
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_sets[i]);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

      // ------------ OR ------------
      info.operation = "OR";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_sets.size(); ++i) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(bit_sets[i]);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(bit_sets[i]);
        
        bs0.Toggle(i%bs0.size());

        bs0 = bs0.OR(bs1);
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (size_t i = 0; i < num_iterations; ++i) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(bit_sets[i]);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;



    }
  }
  // Test bitvectors with std::list as container
  template<size_t NUM_BITS>
  void BenchmarkBitVectors_list() {

    info.container = container_list;

    const size_t width = NUM_BITS;
    info.bits = width;
    // make and resize all bit vectors
    for (size_t i = 0; i < num_iterations; ++i){
      emp::BitVector bitvector1;
      emp::BitVector bitvector2;

      bitvector1.Resize(width);
      bitvector2.Resize(width);

      emp::RandomizeBitVector(bitvector1, random);
      emp::RandomizeBitVector(bitvector2, random);

      auto ntuple = std::make_tuple (bitvector1, bitvector2);


      list_of_bit_vectors.push_back(ntuple);

    }
    // Compute timings for each replicate.
    for (size_t rep = 0; rep < num_replicates; ++rep) {
      info.replicate = rep;
      // Do bit vector:
      info.treatment = "bit_vector";


      emp::BitVector recipient;

      // (2) Time each operation!

      // ------------ CountOnes_Mixed ------------
      info.operation = "CountOnes_Mixed";
      double comulative = 0;

      auto j = 0;
      auto start_time = std::chrono::high_resolution_clock::now();

      for (auto i : list_of_bit_vectors) {
        emp::BitVector & bv = std::get<0>(i);
        bv[j%bv.size()] = !bv[j%bv.size()];
        comulative += bv.CountOnes_Mixed();
        j += 1;
      }

      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();



      // ------------ NOT ------------
      info.operation = "NOT";

      j = 0;
      start_time = std::chrono::high_resolution_clock::now();

      for (auto i : list_of_bit_vectors) {
        emp::BitVector & bv = std::get<0>(i);
        bv[j%bv.size()] = !bv[j%bv.size()];
        bv = bv.NOT();
        j += 1;
      }
      

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
      comulative = 0;

      for (auto i : list_of_bit_vectors) {
        emp::BitVector & bv = std::get<0>(i);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;




      // ------------ AND ------------
      info.operation = "AND";

      j = 0;
      start_time = std::chrono::high_resolution_clock::now();

      for (auto i : list_of_bit_vectors) {
        emp::BitVector & bv0 = std::get<0>(i);
        emp::BitVector & bv1 = std::get<0>(i);
        bv0[j%bv0.size()] = !bv0[j%bv0.size()];
        bv0 = bv0.AND(bv1);
        j += 1;
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      for (auto i : list_of_bit_vectors) {
        emp::BitVector & bv = std::get<0>(i);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;


      // ------------ OR ------------
      info.operation = "OR";

      j = 0;
      start_time = std::chrono::high_resolution_clock::now();
      
      for (auto i : list_of_bit_vectors) {
        emp::BitVector & bv0 = std::get<0>(i);
        emp::BitVector & bv1 = std::get<0>(i);
        bv0[j%bv0.size()] = !bv0[j%bv0.size()];
        bv0 = bv0.OR(bv1);
        j += 1;
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      for (auto i : list_of_bit_vectors) {
        emp::BitVector & bv = std::get<0>(i);
        comulative += bv.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;
    }
  }

  // Test bit sets with list as container
  template<size_t NUM_BITS>
  void BenchmarkBitSets_list(std::list< std::tuple<emp::BitSet<NUM_BITS>,emp::BitSet<NUM_BITS>> > & bit_sets) {
    
    info.container = container_list;
    
    // Resize all bit sets
    const size_t width = NUM_BITS;
    info.bits = width;

    for (size_t i = 0; i < num_iterations; ++i){
      emp::BitSet<width> bitset1;
      emp::BitSet<width> bitset2;


      bitset1.Randomize(random);
      bitset2.Randomize(random);

      auto ntuple = std::make_tuple(bitset1, bitset2);


      bit_sets.push_back(ntuple);

    }

    // Compute timings for each replicate.
    for (size_t rep = 0; rep < num_replicates; ++rep) {
      info.replicate = rep;
      // Do bit set:
      info.treatment = "bit_set";



      double comulative = 0;

      // (2) Time each operation!
      // ------------ CountOnes_Mixed ------------
      info.operation = "CountOnes_Mixed";
      auto start_time = std::chrono::high_resolution_clock::now();
      auto j = 0;
      for (auto i : bit_sets) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(i);
        bs.Toggle(j%bs.size());
        comulative += bs.CountOnes_Mixed();
        j += 1;
      }

      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ NOT ------------
      info.operation = "NOT";
      j = 0;
      start_time = std::chrono::high_resolution_clock::now();

      for (auto i : bit_sets) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(i);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(i);
        
        bs0.Toggle(j%bs0.size());

        bs0 = bs0.NOT();

      j += 1;
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (auto i : bit_sets) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(i);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;


      // ------------ AND ------------
      info.operation = "AND";

      j = 0;
      start_time = std::chrono::high_resolution_clock::now();

      for (auto i : bit_sets) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(i);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(i);
        
        bs0.Toggle(j%bs0.size());

        bs0 = bs0.AND(bs1);

        j += 1;
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      comulative = 0;

      for (auto i : bit_sets) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(i);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;

      // ------------ OR ------------
      info.operation = "OR";
      j = 0;

      start_time = std::chrono::high_resolution_clock::now();

      for (auto i : bit_sets) {
        emp::BitSet<NUM_BITS> & bs0 = std::get<0>(i);
        emp::BitSet<NUM_BITS> & bs1 = std::get<1>(i);
        
        bs0.Toggle(j%bs0.size());

        bs0 = bs0.OR(bs1);

        j += 1;
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      for (auto i : bit_sets) {
        emp::BitSet<NUM_BITS> & bs = std::get<0>(i);
        comulative += bs.CountOnes_Mixed();
      }
      std::cout << comulative << std::endl;
    }

  }
public:
  Benchmark(size_t seed, size_t n_iterations, size_t n_replicates)
    : num_iterations(n_iterations), num_replicates(n_replicates),
      #ifdef EMP_NDEBUG
        datafile("performance.csv"),
      #else
        datafile("performance-debug.csv"),
      #endif
      random(seed),
      bit_vectors(NUM_ITERATIONS),
      bit_set_16(NUM_ITERATIONS),
      bit_set_31(NUM_ITERATIONS),
      bit_set_32(NUM_ITERATIONS),
      bit_set_64(NUM_ITERATIONS),
      bit_set_128(NUM_ITERATIONS),
      bit_set_256(NUM_ITERATIONS),
      bit_set_10000(NUM_ITERATIONS)
  {
    info.n = NUM_ITERATIONS;

    // -- Setup data file --
    // File headers:
    // - rep        => replicate ID
    // - treatment  => treatment ID
    // - operation  => operation
    // - time       => time (s)
    // - n          => number iterations
    datafile.AddFun<size_t>([this]() { return info.replicate; }, "replicate_id");
    datafile.AddFun<std::string>([this]() { return info.treatment; }, "treatment");
    datafile.AddFun<size_t>([this]() { return info.bits; }, "bits");
    datafile.AddFun<std::string>([this]() { return info.operation; }, "operation");
    datafile.AddFun<double>([this]() { return info.time; }, "time");
    datafile.AddFun<size_t>([this]() { return info.n; }, "iterations");
    #ifdef EMP_NDEBUG
      datafile.AddFun<bool>([this]() { return true; }, "compile_time_optimizations");
    #else
      datafile.AddFun<bool>([this]() { return false; }, "compile_time_optimizations");
    #endif
    datafile.AddFun<std::string>([this]() { return info.container; }, "container");

    datafile.PrintHeaderKeys();
  }

  void Run() {
    BenchmarkBitVectors<64>();
    BenchmarkBitSets<64>(lone_bit_set_64);


    BenchmarkBitVectors_vectorized<64>();
    BenchmarkBitSets_vectorized<64>(bit_set_64);
    std::cout << "list of bitvectors" << std::endl;
    BenchmarkBitVectors_list<64>();

    std::cout << "--------------------list of bitsets" << std::endl;

    BenchmarkBitSets_list<64>(list_of_bit_set_64);

    // BenchmarkBitVectors<128>();
    // BenchmarkBitSets<128>(bit_set_128);

    // BenchmarkBitVectors<256>();
    // BenchmarkBitSets<256>(bit_set_256);

    // BenchmarkBitVectors<10000>();
    // BenchmarkBitSets<10000>(bit_set_10000);
  }

};

int main() {
  Benchmark benchmark(SEED, NUM_ITERATIONS, NUM_REPLICATES);
  benchmark.Run();
}