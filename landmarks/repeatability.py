"""Evaluation-only one-to-one candidate proximity, not learned matching."""
import numpy as np
from scipy.optimize import linear_sum_assignment
from scipy.spatial.distance import cdist


def repeatability(a_xy,b_xy,tolerances=(1.,2.,5.)):
    a=np.asarray(a_xy); b=np.asarray(b_xy)
    result=[]
    for tolerance in tolerances:
        if not len(a) or not len(b): count=0
        else:
            distances=cdist(a,b)
            # Large invalid cost makes maximum feasible cardinality primary;
            # proximity breaks ties. No point can be counted twice.
            cost=np.where(distances<=tolerance,distances/max(tolerance,1e-12),(min(len(a),len(b))+1)*2.)
            rows,cols=linear_sum_assignment(cost)
            count=int((distances[rows,cols]<=tolerance).sum())
        result.append(dict(tolerance_mm=float(tolerance),one_to_one_pairs=count,
            fraction_of_smaller_set=count/min(len(a),len(b)) if min(len(a),len(b)) else 0.,
            fraction_of_source=count/len(a) if len(a) else 0.,fraction_of_target=count/len(b) if len(b) else 0.))
    return result
