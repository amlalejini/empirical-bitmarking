#include <iostream>
#include <chrono>
#include <array>
#include <tuple>
#include <map>


#include "base/vector.h"

#include "tools/BitSet.h"
#include "tools/BitVector.h"

#include "data/DataFile.h"

#include "tools/Random.h"
#include "tools/random_utils.h"

#include <thread>

size_t SEED = 2;
size_t NUM_REPLICATES = 100;    // How many replicates should we do?
size_t NUM_ITERATIONS = 100;  // How many times to repeat each operation before measuring time elapsed?


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
  emp::vector< std::tuple<emp::BitSet<64>,emp::BitSet<64>> > bit_set_64;

  const std::string container = "T";
  const std::string container_vectorized = "emp::vector<T>";


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
        std::cout << std::get<0>(bit_vectors[i]) << std::get<1>(bit_vectors[i]) << std::endl;
      }
      emp::BitVector recipient;

      // (2) Time each operation!

      // ------------ CountOnes_Mixed ------------
      info.operation = "CountOnes_Mixed";
      auto start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_vectors.size(); ++i) {
        size_t num_ones = std::get<0>(bit_vectors[i]).CountOnes_Mixed();
      }

      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ NOT ------------
      info.operation = "NOT";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_vectors.size(); ++i) {
         recipient = std::get<0>(bit_vectors[i]).NOT();
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ AND ------------
      info.operation = "AND";
      start_time = std::chrono::high_resolution_clock::now();
   
      for (size_t i = 0; i < bit_vectors.size(); ++i) {
         recipient = std::get<0>(bit_vectors[i]).AND(std::get<1>(bit_vectors[i]));
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ OR ------------
      info.operation = "OR";
      start_time = std::chrono::high_resolution_clock::now();
      
      for (size_t i = 0; i < bit_vectors.size(); ++i) {
         recipient = std::get<0>(bit_vectors[i]).OR(std::get<1>(bit_vectors[i]));
      }
    info.container = container_vectorized;
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
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
        std::cout << std::get<0>(bit_sets[i]) << std::get<1>(bit_sets[i]) << std::endl;
      }
      emp::BitSet<NUM_BITS> recipient;

      // (2) Time each operation!
      // ------------ CountOnes_Mixed ------------
      info.operation = "CountOnes_Mixed";
      auto start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_sets.size(); ++i) {
        size_t num_ones = std::get<0>(bit_sets[i]).CountOnes_Mixed();
      }

      auto end_time = std::chrono::high_resolution_clock::now();
      auto duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ NOT ------------
      info.operation = "NOT";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_sets.size(); ++i) {
        recipient = std::get<0>(bit_sets[i]).NOT();
      }
      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ AND ------------
      info.operation = "AND";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_sets.size(); ++i) {
        recipient = std::get<0>(bit_sets[i]).AND(std::get<1>(bit_sets[i]));

      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();

      // ------------ OR ------------
      info.operation = "OR";
      start_time = std::chrono::high_resolution_clock::now();

      for (size_t i = 0; i < bit_sets.size(); ++i) {
         recipient = std::get<0>(bit_sets[i]).OR(std::get<1>(bit_sets[i]));
      }

      end_time = std::chrono::high_resolution_clock::now();
      duration = std::chrono::duration_cast<std::chrono::nanoseconds>(end_time - start_time).count();
      info.time = duration;
      datafile.Update();
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
      bit_set_64(NUM_ITERATIONS)
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
    BenchmarkBitVectors_vectorized<64>();
    BenchmarkBitSets_vectorized<64>(bit_set_64);

  }
};

int main() {
  Benchmark benchmark(SEED, NUM_ITERATIONS, NUM_REPLICATES);
  benchmark.Run();
}