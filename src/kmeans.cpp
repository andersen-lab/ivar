#include "dataanalysis.h"
#include "ap.h"
#include "stdafx.h"
using namespace alglib;
#include <numeric>
/*
 * Written to extend alglib kmeans clustering, imposing 
 * and additional "compositional constraint" and fixing a noise
 * cluster at 0.03. Additionall changes the scoring metric from 
 * entropy to sil score.
 */
double average(std::vector<double> x){
  double sumTotal = 0;
  for(double k=0; k < x.size(); ++k){
      sumTotal += x[k];
  }
  return(sumTotal / x.size());
}

float calculate_sil_single_point(std::vector<std::vector<double>> sorted_points, std::vector<double> centers, double point, double center){
  /*
  */
  uint32_t n_clusters = centers.size();
  std::vector<double> internal_dist; //keep track of distance between point of interest and all points in cluster
  double external_dist [n_clusters][2]; //dim 1 cluster value, dim 2 distance sum points and point count
  double compare_point=0;
  double a=0; //sil score a term
  double b=1; //sil score b term
  double tmp=0;

  //iterates through every point
  for (uint32_t i = 0; i < sorted_points.size(); i++) {
    for(uint32_t x = 0; x < sorted_points[i].size(); x++){
      compare_point = sorted_points[i][x];
      double cur_center = i;
      //belongs to the same cluster
      if(cur_center == center){
        internal_dist.push_back(abs(point-compare_point));
      }else{ //doesn't belong to the same cluster
        external_dist[i][0] += abs(point-compare_point); //sum distance between poi and cluster points
        external_dist[i][1] += 1; //count number points in cluster
      }
    }
  }
  a = average(internal_dist);
  //find the average distance from each external cluster
  for(uint32_t z=0; z < n_clusters; z++){
    //our target cluster isn't external
    if(z == center){
      continue;
    }
    tmp = external_dist[z][0] / external_dist[z][1];
    //find the minimum of the external cluster dists
    if(tmp < b){
      b = tmp;
    }
  }
  return((b - a) / std::max(a,b));
}
void sil_score(std::vector<std::vector<double>> sorted_points, std::vector<double> centers){
  double index= 0;
  std::vector<double> sil_scores;
  for(std::vector<double> cluster : sorted_points){
    for(double point : cluster){
      sil_scores.push_back(calculate_sil_single_point(sorted_points, centers, point, index));
    }
    index += 1;
  }
  //test lines
  //std::cout << "sil " << average(sil_scores) << std::endl;
}

double find_closest_center(std::vector<double> centers, alglib_impl::ae_int_t k, double point){
  double tmp_distance = 1.0;
  double center;
  double center_index;

  for(int j=0; j<=k-1; j++){
    center = centers[j];
    if(abs(center-point) < tmp_distance){
      tmp_distance = abs(center-point);
      center_index = j;
    }
  }  
  return(center_index);
}

double distance_from_one(std::vector<double> centers){
  double sum = 0;
  double one = 1;
  for(double c : centers){
    sum += c;
  }
  return(one-sum);
}

double calculate_variance(std::vector<double> points){
  double mean = 0.0;
  double accumulator = 0.0;

  mean = std::accumulate(points.begin(), points.end(), 0.0) / points.size();
  for(double point : points){
    accumulator += std::pow(point - mean, 2);
  }
  return(accumulator/points.size());
}

std::vector<double> calculate_variance(std::vector<std::vector<double>> sorted_points){
  std::vector<double> variance;
  for(std::vector<double> points : sorted_points){
    variance.push_back(calculate_variance(points));
  }
  return(variance);
}

std::vector<double> remove_outlier_points(std::vector<double> points, double distance){
  /*
   * Given a cluster, remove the point contributing the most to the variance in the direction
   * opposite the distance and return.
   */ 
  //std::cout << "top of remove outlier, points size : " << points.size() << std::endl;
  std::vector<double> zscores;
  double mean = std::accumulate(points.begin(), points.end(), 0.0) / points.size();
  double stddev = sqrt(calculate_variance(points));

  for(double p: points){
    zscores.push_back((p - mean) / stddev);
  }
  double remove_index = 0;
  if(distance > 0){
    remove_index = std::min_element(zscores.begin(),zscores.end()) - zscores.begin();
      
  }else{
    remove_index = std::max_element(zscores.begin(),zscores.end()) - zscores.begin(); 
  }
  points.erase(points.begin() + remove_index);
  //std::cout << "bottom of remove outlier, points size : " << points.size() << std::endl;
  return(points);

}

std::vector<double> recalculate_centroids(std::vector<std::vector<double>> sorted_points){
  std::vector<double> centers;
  for(std::vector<double> cluster : sorted_points){
    centers.push_back(std::accumulate(cluster.begin(), cluster.end(), 0.0) / cluster.size());
  }
  return(centers);
}

std::vector<double> flatten(std::vector<std::vector<double>> sorted_points){
  std::vector<double> flat_points;
  for(std::vector<double> cluster : sorted_points){
    for(double point : cluster){
      flat_points.push_back(point);
    }
  }
  return(flat_points);
}

std::vector<double> compositional_constraint(alglib_impl::kmeansbuffers *buf, alglib_impl::ae_int_t k, alglib_impl::ae_matrix* xy, alglib_impl::ae_int_t npoints){
  std::vector<double> centers;
  std::vector<double> variance;
  std::vector<double> flat_points;
  std::vector<std::vector<double>> sorted_points(k);
  double point;
  double distance;
  double center_index;
  int max_variance_index;

  for(int j=0; j<=k-1; j++){
    centers.push_back(buf->ct.ptr.pp_double[j][0]);
  }
  //distance from 1
  distance = distance_from_one(centers);
  //initialize values
  for(int i=0; i<=npoints-1; i++){
    point = xy->ptr.pp_double[i][0];
    flat_points.push_back(point);
    center_index = find_closest_center(centers, k, point);
    sorted_points[center_index].push_back(point);
  }

  uint32_t minimum_flat_points = flat_points.size() * 0.25;
  bool done = false;

  while((abs(distance) > 0.03) && (!done) && (flat_points.size() > minimum_flat_points)){
  //while(abs(distance) > 0.03){
   //calculate the variance per cluster
    variance = calculate_variance(sorted_points);
    max_variance_index = std::max_element(variance.begin(),variance.end()) - variance.begin();

    //remove the point contributing most to variance in the wrong direction
    sorted_points[max_variance_index] = remove_outlier_points(sorted_points[max_variance_index], distance);

    //recalculate centroid
    centers.clear();
    centers = recalculate_centroids(sorted_points);

    //flatten points out again
    flat_points.clear();
    flat_points = flatten(sorted_points);

    //recalculate distance
    distance = distance_from_one(centers); 
    
    //try and determine if we have any clusters with few points or if we've lost too many points total
    for(std::vector<double> cluster : sorted_points){
      if(cluster.size() < 2){
        done = true;
      }
    }

  }
  return(centers);  
}

static ae_bool clustering_fixcenters(alglib_impl::ae_matrix* xy,
     alglib_impl::ae_int_t npoints,
     alglib_impl::ae_int_t nvars,
     alglib_impl::ae_matrix* ct,
     alglib_impl::ae_int_t k,
     alglib_impl::apbuffers* initbuf,
     alglib_impl::ae_shared_pool* updatepool,
     alglib_impl::ae_state *_state)
{
    alglib_impl::ae_int_t fixiteration;
    alglib_impl::ae_int_t centertofix;
    alglib_impl::ae_int_t i;
    alglib_impl::ae_int_t j;
    alglib_impl::ae_int_t pdistant;
    double ddistant;
    double v;
    ae_bool result;


    alglib_impl::ae_assert(npoints>=1, "FixCenters: internal error", _state);
    alglib_impl::ae_assert(nvars>=1, "FixCenters: internal error", _state);
    alglib_impl::ae_assert(k>=1, "FixCenters: internal error", _state);

    /*
     * Calculate distances from points to best centers (RA0)
     * and best center indexes (IA0)
     */
    ivectorsetlengthatleast(&initbuf->ia0, npoints, _state);
    rvectorsetlengthatleast(&initbuf->ra0, npoints, _state);
    alglib_impl::kmeansupdatedistances(xy, 0, npoints, nvars, ct, 0, k, &initbuf->ia0, &initbuf->ra0, updatepool, _state);

    /*
     * Repeat loop:
     * * find first center which has no corresponding point
     * * set it to the most distant (from the rest of the centerset) point
     * * recalculate distances, update IA0/RA0
     * * repeat
     *
     * Loop is repeated for at most 2*K iterations. It is stopped once we have
     * no "empty" clusters.
     */
    bvectorsetlengthatleast(&initbuf->ba0, k, _state);
    for(fixiteration=0; fixiteration<=2*k; fixiteration++)
    {
        /*
         * Select center to fix (one which is not mentioned in IA0),
         * terminate if there is no such center.
         * BA0[] stores True for centers which have at least one point.
         */
        for(i=0; i<=k-1; i++)
        {
            initbuf->ba0.ptr.p_bool[i] = ae_false;
        }
        for(i=0; i<=npoints-1; i++)
        {
            initbuf->ba0.ptr.p_bool[initbuf->ia0.ptr.p_int[i]] = ae_true;
        }
        centertofix = -1;
        for(i=0; i<=k-1; i++)
        {
            if( !initbuf->ba0.ptr.p_bool[i] )
            {
                centertofix = i;
                break;
            }
        }
        if( centertofix<0 )
        {
            result = ae_true;
            return result;
        }

        /*
         * Replace center to fix by the most distant point.
         * Update IA0/RA0
         */
        pdistant = 0;
        ddistant = initbuf->ra0.ptr.p_double[pdistant];
        for(i=0; i<=npoints-1; i++)
        {
            if( alglib_impl::ae_fp_greater(initbuf->ra0.ptr.p_double[i],ddistant) )
            {
                ddistant = initbuf->ra0.ptr.p_double[i];
                pdistant = i;
            }
        }
        if( alglib_impl::ae_fp_eq(ddistant,0.0) )
        {
            break;
        }
        alglib_impl::ae_v_move(&ct->ptr.pp_double[centertofix][0], 1, &xy->ptr.pp_double[pdistant][0], 1, alglib_impl::ae_v_len(0,nvars-1));
        for(i=0; i<=npoints-1; i++)
        {
            v = 0.0;
            for(j=0; j<=nvars-1; j++)
            {
                v = v+ae_sqr(xy->ptr.pp_double[i][j]-ct->ptr.pp_double[centertofix][j], _state);
            }
            if( alglib_impl::ae_fp_less(v,initbuf->ra0.ptr.p_double[i]) )
            {
                initbuf->ra0.ptr.p_double[i] = v;
                initbuf->ia0.ptr.p_int[i] = centertofix;
            }
        }
    }
    result = ae_false;
    return result;
}
static void clustering_selectinitialcenters(alglib_impl::ae_matrix* xy,
     alglib_impl::ae_int_t npoints,
     alglib_impl::ae_int_t nvars,
     alglib_impl::ae_int_t initalgo,
     alglib_impl::hqrndstate* rs,
     alglib_impl::ae_int_t k,
     alglib_impl::ae_matrix* ct,
     alglib_impl::apbuffers* initbuf,
     alglib_impl::ae_shared_pool* updatepool,
     alglib_impl::ae_state *_state)
{
    alglib_impl::ae_int_t cidx;
    alglib_impl::ae_int_t i;
    alglib_impl::ae_int_t j;
    double v;
    double vv;
    double s;
    alglib_impl::ae_int_t lastnz;
    alglib_impl::ae_int_t ptidx;
    alglib_impl::ae_int_t samplesize;
    alglib_impl::ae_int_t samplescntnew;
    alglib_impl::ae_int_t samplescntall;
    double samplescale;

    /*
     * Check parameters
     */
    alglib_impl::ae_assert(npoints>0, "SelectInitialCenters: internal error", _state);
    alglib_impl::ae_assert(nvars>0, "SelectInitialCenters: internal error", _state);
    alglib_impl::ae_assert(k>0, "SelectInitialCenters: internal error", _state);
    if( initalgo==0 )
    {
        initalgo = 3;
    }
    rmatrixsetlengthatleast(ct, k, nvars, _state);

    /*
     * Random initialization
     */
    if( initalgo==-1 )
    {
        for(i=0; i<=k-1; i++)
        {
          alglib_impl::ae_v_move(&ct->ptr.pp_double[i][0], 1, &xy->ptr.pp_double[i%npoints][0], 1, alglib_impl::ae_v_len(0,nvars-1));
        }
        return;
    }

    /*
     * Random initialization
     */
        if( initalgo==1 )
    {
        for(i=0; i<=k-1; i++)
        {
            j = hqrnduniformi(rs, npoints, _state);
            alglib_impl::ae_v_move(&ct->ptr.pp_double[i][0], 1, &xy->ptr.pp_double[j][0], 1, alglib_impl::ae_v_len(0,nvars-1));
        }
        return;
    }

    /*
     * k-means++ initialization
     */
    if( initalgo==2 )
    {

        /*
         * Prepare distances array.
         * Select initial center at random.
         */
        rvectorsetlengthatleast(&initbuf->ra0, npoints, _state);
        for(i=0; i<=npoints-1; i++)
        {
            initbuf->ra0.ptr.p_double[i] = ae_maxrealnumber;
        }
        ptidx = hqrnduniformi(rs, npoints, _state);
        alglib_impl::ae_v_move(&ct->ptr.pp_double[0][0], 1, &xy->ptr.pp_double[ptidx][0], 1, alglib_impl::ae_v_len(0,nvars-1));

        /*
         * For each newly added center repeat:
         * * reevaluate distances from points to best centers
         * * sample points with probability dependent on distance
         * * add new center
         */
        for(cidx=0; cidx<=k-2; cidx++)
        {

            /*
             * Reevaluate distances
             */
            s = 0.0;
            for(i=0; i<=npoints-1; i++)
            {
                v = 0.0;
                for(j=0; j<=nvars-1; j++)
                {
                    vv = xy->ptr.pp_double[i][j]-ct->ptr.pp_double[cidx][j];
                    v = v+vv*vv;
                }
                if( alglib_impl::ae_fp_less(v,initbuf->ra0.ptr.p_double[i]) )
                {
                    initbuf->ra0.ptr.p_double[i] = v;
                }
                s = s+initbuf->ra0.ptr.p_double[i];
            }
                        /*
             * If all distances are zero, it means that we can not find enough
             * distinct points. In this case we just select non-distinct center
             * at random and continue iterations. This issue will be handled
             * later in the FixCenters() function.
             */
            if( alglib_impl::ae_fp_eq(s,0.0) )
            {
                ptidx = alglib_impl::hqrnduniformi(rs, npoints, _state);
                alglib_impl::ae_v_move(&ct->ptr.pp_double[cidx+1][0], 1, &xy->ptr.pp_double[ptidx][0], 1, alglib_impl::ae_v_len(0,nvars-1));
                continue;
            }

            /*
             * Select point as center using its distance.
             * We also handle situation when because of rounding errors
             * no point was selected - in this case, last non-zero one
             * will be used.
             */
            v = alglib_impl::hqrnduniformr(rs, _state);
            vv = 0.0;
            lastnz = -1;
            ptidx = -1;
            for(i=0; i<=npoints-1; i++)
            {
                if(alglib_impl::ae_fp_eq(initbuf->ra0.ptr.p_double[i],0.0) )
                {
                    continue;
                }
                lastnz = i;
                vv = vv+initbuf->ra0.ptr.p_double[i];
                if(alglib_impl::ae_fp_less_eq(v,vv/s) )
                {
                    ptidx = i;
                    break;
                }
            }
            alglib_impl::ae_assert(lastnz>=0, "SelectInitialCenters: integrity error", _state);
            if( ptidx<0 )
            {
                ptidx = lastnz;
            }
            alglib_impl::ae_v_move(&ct->ptr.pp_double[cidx+1][0], 1, &xy->ptr.pp_double[ptidx][0], 1, alglib_impl::ae_v_len(0,nvars-1));
        }
        return;
    }
    /*
     * "Fast-greedy" algorithm based on "Scalable k-means++".
     *
     * We perform several rounds, within each round we sample about 0.5*K points
     * (not exactly 0.5*K) until we have 2*K points sampled. Before each round
     * we calculate distances from dataset points to closest points sampled so far.
     * We sample dataset points independently using distance xtimes 0.5*K divided by total
     * as probability (similar to k-means++, but each point is sampled independently;
     * after each round we have roughtly 0.5*K points added to sample).
     *
     * After sampling is done, we run "greedy" version of k-means++ on this subsample
     * which selects most distant point on every round.
     */
    if( initalgo==3 )
    {

        /*
         * Prepare arrays.
         * Select initial center at random, add it to "new" part of sample,
         * which is stored at the beginning of the array
         */
        samplesize = 2*k;
        samplescale = 0.5*k;
        rmatrixsetlengthatleast(&initbuf->rm0, samplesize, nvars, _state);
        ptidx = hqrnduniformi(rs, npoints, _state);
        alglib_impl::ae_v_move(&initbuf->rm0.ptr.pp_double[0][0], 1, &xy->ptr.pp_double[ptidx][0], 1, alglib_impl::ae_v_len(0,nvars-1));
        samplescntnew = 1;
        samplescntall = 1;
        rvectorsetlengthatleast(&initbuf->ra0, npoints, _state);
        rvectorsetlengthatleast(&initbuf->ra1, npoints, _state);
        ivectorsetlengthatleast(&initbuf->ia1, npoints, _state);
        for(i=0; i<=npoints-1; i++)
        {
            initbuf->ra0.ptr.p_double[i] = ae_maxrealnumber;
        }

        /*
         * Repeat until samples count is 2*K
         */
        while(samplescntall<samplesize)
        {

            /*
             * Evaluate distances from points to NEW centers, store to RA1.
             * Reset counter of "new" centers.
             */
            kmeansupdatedistances(xy, 0, npoints, nvars, &initbuf->rm0, samplescntall-samplescntnew, samplescntall, &initbuf->ia1, &initbuf->ra1, updatepool, _state);
            samplescntnew = 0;
                        /*
             * Merge new distances with old ones.
             * Calculate sum of distances, if sum is exactly zero - fill sample
             * by randomly selected points and terminate.
             */
            s = 0.0;
            for(i=0; i<=npoints-1; i++)
            {
                initbuf->ra0.ptr.p_double[i] = ae_minreal(initbuf->ra0.ptr.p_double[i], initbuf->ra1.ptr.p_double[i], _state);
                s = s+initbuf->ra0.ptr.p_double[i];
            }
            if( alglib_impl::ae_fp_eq(s,0.0) )
            {
                while(samplescntall<samplesize)
                {
                    ptidx = hqrnduniformi(rs, npoints, _state);
                    alglib_impl::ae_v_move(&initbuf->rm0.ptr.pp_double[samplescntall][0], 1, &xy->ptr.pp_double[ptidx][0], 1, alglib_impl::ae_v_len(0,nvars-1));
                    inc(&samplescntall, _state);
                    inc(&samplescntnew, _state);
                }
                break;
            }

            /*
             * Sample points independently.
             */
            for(i=0; i<=npoints-1; i++)
            {
                if( samplescntall==samplesize )
                {
                    break;
                }
                if( alglib_impl::ae_fp_eq(initbuf->ra0.ptr.p_double[i],0.0) )
                {
                    continue;
                }
                if( alglib_impl::ae_fp_less_eq(hqrnduniformr(rs, _state),samplescale*initbuf->ra0.ptr.p_double[i]/s) )
                {
                  alglib_impl::ae_v_move(&initbuf->rm0.ptr.pp_double[samplescntall][0], 1, &xy->ptr.pp_double[i][0], 1, alglib_impl::ae_v_len(0,nvars-1));
                    inc(&samplescntall, _state);
                    inc(&samplescntnew, _state);
                }
            }
        }

        /*
         * Run greedy version of k-means on sampled points
         */
                rvectorsetlengthatleast(&initbuf->ra0, samplescntall, _state);
        for(i=0; i<=samplescntall-1; i++)
        {
            initbuf->ra0.ptr.p_double[i] = ae_maxrealnumber;
        }
        ptidx = hqrnduniformi(rs, samplescntall, _state);
        alglib_impl::ae_v_move(&ct->ptr.pp_double[0][0], 1, &initbuf->rm0.ptr.pp_double[ptidx][0], 1, alglib_impl::ae_v_len(0,nvars-1));
        for(cidx=0; cidx<=k-2; cidx++)
        {

            /*
             * Reevaluate distances
             */
            for(i=0; i<=samplescntall-1; i++)
            {
                v = 0.0;
                for(j=0; j<=nvars-1; j++)
                {
                    vv = initbuf->rm0.ptr.pp_double[i][j]-ct->ptr.pp_double[cidx][j];
                    v = v+vv*vv;
                }
                if( alglib_impl::ae_fp_less(v,initbuf->ra0.ptr.p_double[i]) )
                {
                    initbuf->ra0.ptr.p_double[i] = v;
                }
            }

            /*
             * Select point as center in greedy manner - most distant
             * point is selected.
             */
            ptidx = 0;
            for(i=0; i<=samplescntall-1; i++)
            {
                if( alglib_impl::ae_fp_greater(initbuf->ra0.ptr.p_double[i],initbuf->ra0.ptr.p_double[ptidx]) )
                {
                    ptidx = i;
                }
            }
            alglib_impl::ae_v_move(&ct->ptr.pp_double[cidx+1][0], 1, &initbuf->rm0.ptr.pp_double[ptidx][0], 1, alglib_impl::ae_v_len(0,nvars-1));
        }
        return;
    }
    /*
     * Internal error
     */
    alglib_impl::ae_assert(ae_false, "SelectInitialCenters: internal error", _state);
}

void kmeans_internal(alglib_impl::ae_matrix* xy, alglib_impl::ae_int_t npoints, alglib_impl::ae_int_t nvars, alglib_impl::ae_int_t k, alglib_impl::ae_int_t initalgo, alglib_impl::ae_int_t seed, alglib_impl::ae_int_t maxits, alglib_impl::ae_int_t restarts, ae_bool kmeansdbgnoits, alglib_impl::ae_int_t* info, alglib_impl::ae_int_t* iterationscount, alglib_impl::ae_matrix* ccol, ae_bool needccol, alglib_impl::ae_matrix* crow, ae_bool needcrow, alglib_impl::ae_vector* xyc, double* energy, alglib_impl::kmeansbuffers* buf, alglib_impl::ae_state *_state)
{
    alglib_impl::ae_frame _frame_block;
    alglib_impl::ae_int_t i;
    alglib_impl::ae_int_t j;
    alglib_impl::ae_int_t i1;
    double e;
    double eprev;
    double v;
    double vv;
    ae_bool waschanges;
    ae_bool zerosizeclusters;
    alglib_impl::ae_int_t pass;
    alglib_impl::ae_int_t itcnt;
    alglib_impl::hqrndstate rs;

    alglib_impl::ae_frame_make(_state, &_frame_block);
    memset(&rs, 0, sizeof(rs));
    *info = 0;
    *iterationscount = 0;
    alglib_impl::ae_matrix_clear(ccol);
    alglib_impl::ae_matrix_clear(crow);
    alglib_impl::ae_vector_clear(xyc);
    *energy = 0;
    alglib_impl::_hqrndstate_init(&rs, _state, ae_true);
    /*
     * Test parameters
     */
    if( ((npoints<k||nvars<1)||k<1)||restarts<1 )
    {
        *info = -1;
        *iterationscount = 0;
        alglib_impl::ae_frame_leave(_state);
        return;
    }
        /*
     * TODO: special case K=1
     * TODO: special case K=NPoints
     */
    *info = 1;
    *iterationscount = 0;

    /*
     * Multiple passes of k-means++ algorithm
     */
    //make sure the seed is set
    if( seed<=0 )
    {
      alglib_impl::hqrndrandomize(&rs, _state);
    }
    else
    {
      alglib_impl::hqrndseed(325355, seed, &rs, _state);
    }
    alglib_impl::ae_vector_set_length(xyc, npoints, _state);
    rmatrixsetlengthatleast(&buf->ct, k, nvars, _state);
    rmatrixsetlengthatleast(&buf->ctbest, k, nvars, _state);
    ivectorsetlengthatleast(&buf->xycprev, npoints, _state);
    ivectorsetlengthatleast(&buf->xycbest, npoints, _state);
    rvectorsetlengthatleast(&buf->d2, npoints, _state);
    ivectorsetlengthatleast(&buf->csizes, k, _state);
    *energy = ae_maxrealnumber;
    for(pass=1; pass<=restarts; pass++)
    {
        /*
         * Select initial centers.
         *
         * Note that for performance reasons centers are stored in ROWS of CT, not
         * in columns. We'll transpose CT in the end and store it in the C.
         *
         * Also note that SelectInitialCenters() may return degenerate set of centers
         * (some of them have no corresponding points in dataset, some are non-distinct).
         * Algorithm below is robust enough to deal with such set.
         */

        //CHANGE: select the centers from the list of points, fix a cluster at 0.03
        clustering_selectinitialcenters(xy, npoints, nvars, initalgo, &rs, k, &buf->ct, &buf->initbuf, &buf->updatepool, _state);
        /*
         * Lloyd's iteration
         */
        //boolean value, not sure what this means but this executes
        if( !kmeansdbgnoits )
        {
            /*
             * Perform iteration as usual, in normal mode
             */
            for(i=0; i<=npoints-1; i++)
            {
                xyc->ptr.p_int[i] = -1;
            }
            eprev = ae_maxrealnumber;
            e = ae_maxrealnumber;
            itcnt = 0;
            //general notes
            //xyc is are the center values
            while(maxits==0||itcnt<maxits)
            {
                /*
                 * Update iteration counter
                 */
                itcnt = itcnt+1;
                inc(iterationscount, _state);

                /*
                 * Call KMeansUpdateDistances(), fill XYC with center numbers,
                 * D2 with center distances.
                 */
                for(i=0; i<=npoints-1; i++)
                {
                    buf->xycprev.ptr.p_int[i] = xyc->ptr.p_int[i];
                }
                alglib_impl::kmeansupdatedistances(xy, 0, npoints, nvars, &buf->ct, 0, k, xyc, &buf->d2, &buf->updatepool, _state);
                waschanges = ae_false;
                for(i=0; i<=npoints-1; i++)
                {
                    waschanges = waschanges||xyc->ptr.p_int[i]!=buf->xycprev.ptr.p_int[i];
                }

                /*
                 * Update centers
                 */

                //CHANGE: at this point we try and force the centers to equal 1
                for(j=0; j<=k-1; j++)
                {
                    buf->csizes.ptr.p_int[j] = 0;
                }
                for(i=0; i<=k-1; i++)
                {
                    for(j=0; j<=nvars-1; j++)
                    {
                        buf->ct.ptr.pp_double[i][j] = (double)(0);
                    }
                }
                for(i=0; i<=npoints-1; i++)
                {
                    buf->csizes.ptr.p_int[xyc->ptr.p_int[i]] = buf->csizes.ptr.p_int[xyc->ptr.p_int[i]]+1;
                    //vdst, stride_dist, vscr, stride_src, n
                    alglib_impl::ae_v_add(&buf->ct.ptr.pp_double[xyc->ptr.p_int[i]][0], 1, &xy->ptr.pp_double[i][0], 1, alglib_impl::ae_v_len(0,nvars-1));
                }
                                
                //check for orphan clusters
                zerosizeclusters = ae_false;
                for(j=0; j<=k-1; j++)
                {
                    if( buf->csizes.ptr.p_int[j]!=0 )
                    {
                        v = (double)1/(double)buf->csizes.ptr.p_int[j];
                        //this is the line that sets the new values
                        alglib_impl::ae_v_muld(&buf->ct.ptr.pp_double[j][0], 1, alglib_impl::ae_v_len(0,nvars-1), v);
                   }
                    zerosizeclusters = zerosizeclusters||buf->csizes.ptr.p_int[j]==0;
                }
                //enforce compositional constraint
                std::vector<double> centers;
                std::vector<std::vector<double>> sorted_points(k);
                centers = compositional_constraint(buf, k, xy, npoints);
                for(j=0; j<=npoints-1;j++){
                  float c = find_closest_center(centers, k, xy->ptr.pp_double[j][0]);
                  sorted_points[c].push_back(xy->ptr.pp_double[j][0]);
                }
                sil_score(sorted_points, centers);
                
                //fix noise cluster
                int min_index  = std::min_element(centers.begin(), centers.end()) - centers.begin();
                centers[min_index] = 0.03;

                for(j=0; j<=k-1; j++){
                  buf->ct.ptr.pp_double[j][0] = centers[j];
                }

                if( zerosizeclusters )
                {

                    /*
                     * Some clusters have zero size - rare, but possible.
                     * We'll choose new centers for such clusters using k-means++ rule
                     * and restart algorithm, decrementing iteration counter
                     * in order to allow one more iteration (this one was useless
                     * and should not be counted).
                     */
                    if( !clustering_fixcenters(xy, npoints, nvars, &buf->ct, k, &buf->initbuf, &buf->updatepool, _state) )
                    {
                        *info = -3;
                        alglib_impl::ae_frame_leave(_state);
                        return;
                    }
                    itcnt = itcnt-1;
                    continue;
                }
                /*
                 * Stop if one of two conditions is met:
                 * 1. nothing has changed during iteration
                * 2. energy function increased after recalculation on new centers
                */
                e = (double)(0);
                 
                //here we replace energy with sil score
                for(i=0; i<=npoints-1; i++)
                {
                    v = 0.0;
                    i1 = xyc->ptr.p_int[i];
                    for(j=0; j<=nvars-1; j++)
                    {
                        vv = xy->ptr.pp_double[i][j]-buf->ct.ptr.pp_double[i1][j];
                        v = v+vv*vv;
                    }
                    e = e+v;
                }
                if( !waschanges|| alglib_impl::ae_fp_greater_eq(e,eprev) )
                {
                    break;
                }

                /*
                 * Update EPrev
                 */
                eprev = e;
            }
        }
        else
        {

            /*
             * Debug mode: no Lloyd's iteration.
             * We just calculate potential E.
             */
          alglib_impl::kmeansupdatedistances(xy, 0, npoints, nvars, &buf->ct, 0, k, xyc, &buf->d2, &buf->updatepool, _state);
            e = (double)(0);
            for(i=0; i<=npoints-1; i++)
            {
                e = e+buf->d2.ptr.p_double[i];
            }
        }

        /*
         * Compare E with best centers found so far
         */
        if( alglib_impl::ae_fp_less(e,*energy) )
        {

            /*
             * store partition.
             */
            *energy = e;
            copymatrix(&buf->ct, 0, k-1, 0, nvars-1, &buf->ctbest, 0, k-1, 0, nvars-1, _state);
            for(i=0; i<=npoints-1; i++)
            {
                buf->xycbest.ptr.p_int[i] = xyc->ptr.p_int[i];
            }
        }
    }
        /*
     * Copy and transpose
     */
    if( needccol )
    {
      alglib_impl::ae_matrix_set_length(ccol, nvars, k, _state);
        copyandtranspose(&buf->ctbest, 0, k-1, 0, nvars-1, ccol, 0, nvars-1, 0, k-1, _state);
    }
    if( needcrow )
    {
      alglib_impl::ae_matrix_set_length(crow, k, nvars, _state);
        rmatrixcopy(k, nvars, &buf->ctbest, 0, 0, crow, 0, 0, _state);
    }
    for(i=0; i<=npoints-1; i++)
    {
        xyc->ptr.p_int[i] = buf->xycbest.ptr.p_int[i];
    }
    alglib_impl::ae_frame_leave(_state);
}

void custom_kmeans_inner(alglib_impl::clusterizerstate* s, alglib_impl::ae_int_t k, alglib_impl::kmeansreport* rep, alglib_impl::ae_state *_state)
{
    alglib_impl::ae_frame _frame_block;
    alglib_impl::ae_matrix dummy;

    alglib_impl::ae_frame_make(_state, &_frame_block);
    memset(&dummy, 0, sizeof(dummy));
    alglib_impl::_kmeansreport_clear(rep);
    alglib_impl::ae_matrix_init(&dummy, 0, 0, alglib_impl::DT_REAL, _state, ae_true);

    alglib_impl::ae_assert(k>=0, "ClusterizerRunKMeans: K<0", _state);

    /*
     * Incorrect distance type
     */
    if( s->disttype!=2 )
    {
        rep->npoints = s->npoints;
        rep->terminationtype = -5;
        rep->k = k;
        rep->iterationscount = 0;
        rep->energy = 0.0;
        alglib_impl::ae_frame_leave(_state);
        return;
    }

    /*
     * K>NPoints or (K=0 and NPoints>0)
     */
    if( k>s->npoints||(k==0&&s->npoints>0) )
    {
        rep->npoints = s->npoints;
        rep->terminationtype = -3;
        rep->k = k;
        rep->iterationscount = 0;
        rep->energy = 0.0;
        alglib_impl::ae_frame_leave(_state);
        return;
    }

    /*
     * No points
     */
    if( s->npoints==0 )
    {
        rep->npoints = 0;
        rep->terminationtype = 1;
        rep->k = k;
        rep->iterationscount = 0;
        rep->energy = 0.0;
        alglib_impl::ae_frame_leave(_state);
        return;
    }
        /*
     * Normal case:
     * 1<=K<=NPoints, Euclidean distance
     */
    rep->npoints = s->npoints;
    rep->nfeatures = s->nfeatures;
    rep->k = k;
    rep->npoints = s->npoints;
    rep->nfeatures = s->nfeatures;
    //call to function that actually does the clustering 
    kmeans_internal(&s->xy, s->npoints, s->nfeatures, k, s->kmeansinitalgo, s->seed, s->kmeansmaxits, s->kmeansrestarts, s->kmeansdbgnoits, &rep->terminationtype, &rep->iterationscount, &dummy, ae_false, &rep->c, ae_true, &rep->cidx, &rep->energy, &s->kmeanstmp, _state);
    alglib_impl::ae_frame_leave(_state);
}
void custom_kmeans(const alglib::clusterizerstate &s, const alglib_impl::ae_int_t k, alglib::kmeansreport &rep, const xparams _xparams)
{
    jmp_buf _break_jump;
    alglib_impl::ae_state _alglib_env_state;
    alglib_impl::ae_state_init(&_alglib_env_state);
    if( setjmp(_break_jump) )
    {
#if !defined(AE_NO_EXCEPTIONS)
        _ALGLIB_CPP_EXCEPTION(_alglib_env_state.error_msg);
#else
        _ALGLIB_SET_ERROR_FLAG(_alglib_env_state.error_msg);
        return;
#endif
    }
    ae_state_set_break_jump(&_alglib_env_state, &_break_jump);
    if( _xparams.flags!=0x0 )
        ae_state_set_flags(&_alglib_env_state, _xparams.flags);
    custom_kmeans_inner(const_cast<alglib_impl::clusterizerstate*>(s.c_ptr()), k, const_cast<alglib_impl::kmeansreport*>(rep.c_ptr()), &_alglib_env_state);
    alglib_impl::ae_state_clear(&_alglib_env_state);
    return;
}


