# Social Round-Robin Adapter

This folder is a separate sociology-facing adapter for applying Cell-GPS-style
searcher/findee topology to social round-robin rating data.

It is intentionally kept outside `src/cellgps` and `src/sfplot`.

## Mapping

| Round-robin / SRM term | Cell-GPS term |
| --- | --- |
| actor / perceiver / rater | searcher |
| partner / target / rated person | findee |
| score from actor i to partner j | searcher-findee affinity |

Most social scores are affinities: a larger value usually means more liking,
trust, attraction, similarity, or positive evaluation. Cell-GPS distance
matrices use the opposite convention, where smaller values mean closer. This
adapter therefore keeps two matrices:

- `affinity_matrix`: original directed score matrix, rows are searchers and
  columns are findees.
- `distance_matrix`: Cell-GPS-ready directed distance matrix.

For high-is-close ratings, the normalized conversion is:

```text
distance_ij = 1 - normalized(score_ij)
```

The round-robin diagonal is usually missing because people do not rate
themselves. By default this adapter fills the diagonal with a neutral distance
instead of zero, avoiding a fake self-distance signal.

## Minimal Use

```python
import pandas as pd

from social_round_robin import compute_round_robin_cellgps

ratings = pd.DataFrame(
    {
        "actor": ["A", "A", "B", "B", "C", "C"],
        "partner": ["B", "C", "A", "C", "A", "B"],
        "liking": [5, 1, 4, 2, 1, 5],
    }
)

result = compute_round_robin_cellgps(
    ratings,
    actor_col="actor",
    partner_col="partner",
    score_col="liking",
    min_score=1,
    max_score=5,
)

print(result.affinity_matrix)
print(result.distance_matrix)
print(result.row_cophenetic)  # searcher/actor pattern topology
print(result.col_cophenetic)  # findee/partner pattern topology
```

Run from this repository with:

```powershell
$env:PYTHONPATH = "src"
python -m social_round_robin.examples.toy_round_robin
```

## Interpretation

- `row_cophenetic`: clusters people by how similarly they rate others.
- `col_cophenetic`: clusters people by how similarly they are rated by others.

For SRM-style work, a useful next step is to run the adapter on multiple
matrices:

- raw rating matrix
- actor-effect adjusted matrix
- partner-effect adjusted matrix
- relationship residual matrix

The residual matrix is often the most interesting for dyad-specific topology,
because it removes general rater severity and general target popularity.
