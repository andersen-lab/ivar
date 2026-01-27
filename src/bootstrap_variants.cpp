#include "bootstrap_variants.h"

bootstrap_variants::bootstrap_variants(const std::vector<uint32_t> &sites_, const std::vector<uint32_t> &depths_, const std::vector<uint32_t> &total_depths_, unsigned int seed)
 :seed(seed), rng(seed), sites(sites_), depths(depths_), total_depths(total_depths_) {

  if (sites_.size() != depths_.size() || sites_.size() != total_depths_.size()) {
    throw std::invalid_argument("sites, depths, and total_depths must have the same size");
  }

  for (size_t i = 0; i < sites_.size(); ++i) {
    if (depths_[i] > total_depths_[i]) {
      throw std::invalid_argument("depth cannot exceed total_depth at position " + std::to_string(sites_[i]));
    }
  }

  std::unordered_map<uint32_t, uint32_t> site_total_depth_check;
  for (size_t i = 0; i < sites_.size(); ++i) {
    if (site_total_depth_check.count(sites_[i]) > 0) {
      if (site_total_depth_check[sites_[i]] != total_depths_[i]) {
        throw std::invalid_argument("All variants at position " + std::to_string(sites_[i]) + " must have the same total_depth");
      }
    } else {
      site_total_depth_check[sites_[i]] = total_depths_[i];
    }
  }

  this->build_site_index();
  // Setup site_sampler to sample from multinomial weighted by total depth at each site
  this->site_sampler = new multinomial_sampler(this->unique_sites_total_depths, this->seed);
}

void bootstrap_variants::build_site_index() {
  site_index.clear();

  for (int i = 0; i < this->sites.size(); ++i) {
    site_index[this->sites[i]].push_back(i);
  }

  // Extract unique positions and their weights
  this->unique_sites.reserve(this->site_index.size());
  this->unique_sites_total_depths.reserve(this->site_index.size());

  for (const std::pair<const uint32_t, std::vector<uint32_t>>& elem : site_index) {
    this->unique_sites.push_back(elem.first);
    // All variants at same position have same TOTAL_DEPTH
    this->unique_sites_total_depths.push_back(total_depths[elem.second[0]]);
  }
}

void bootstrap_variants::sample(std::vector<uint32_t> &sampled_variants, std::vector<uint32_t> &sampled_depths, std::vector<double> &sampled_frequencies, uint32_t n_sites) {
  sampled_variants.clear();
  sampled_depths.clear();
  sampled_frequencies.clear();

  std::vector<uint32_t> sampled_site_indices;

  // Sample site indices
  sample_sites(sampled_site_indices, n_sites);

  uint32_t new_site_count = 0;
  for(int s = 0; s < sampled_site_indices.size(); s++) {
    uint32_t sampled_site_idx = sampled_site_indices[s];
    uint32_t sampled_site = this->unique_sites[sampled_site_idx];
    std::vector<uint32_t> variant_indices = site_index[sampled_site];
    std::vector<uint32_t> variant_depths;

    variant_depths.reserve(variant_indices.size());

    for (uint32_t var_index : variant_indices) {
      variant_depths.push_back(this->depths[var_index]);
    }

    // Set up variant_sampler to sample from multinomial weighted by variant depths
    multinomial_sampler variant_sampler(variant_depths, this->seed + sampled_site_idx + s);

    uint32_t site_total_depth = this->total_depths[variant_indices[0]];

    // Sample variant indices
    std::vector<uint32_t> sampled_variant_indices;
    variant_sampler.sample(sampled_variant_indices, site_total_depth);

    // Count total depth of each variant_idx
    std::vector<int> sampled_site_variant_depths(variant_depths.size(), 0);
    for (int x : sampled_variant_indices)
      sampled_site_variant_depths[x]++;

    for (uint32_t i = 0; i < sampled_site_variant_depths.size(); ++i) {
      if (sampled_site_variant_depths[i] == 0)
        continue; // skip variants with zero sampeld depth
//      sampled_variants.push_back(variant_indices[i]); // index of sampled variant
      sampled_variants.push_back(new_site_count); // create new site index
      sampled_depths.push_back(sampled_site_variant_depths[i]); // sampled depth of variant
      double freq = static_cast<double>(sampled_site_variant_depths[i]) / static_cast<double>(site_total_depth);
      sampled_frequencies.push_back(freq);
      new_site_count++;
    }
  }
}

void bootstrap_variants::sample_sites(std::vector<uint32_t> &sampled_sites, uint32_t n_sites) {
  this->site_sampler->sample(sampled_sites, n_sites);
}