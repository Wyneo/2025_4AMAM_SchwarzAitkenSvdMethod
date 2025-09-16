# Décomposition de domaines de Schwarz accélérée par la méthode Aitken-SVD

### Mots-clés : Analyse numérique, Calcul haute performance, Décomposition de domaines, Méthodes d'accélération de convergence, Équations aux dérivées partielles, Méthode d'éléments finis de Galerkine, Matlab

### Projet : Vérifier l'accélération d'Aitken-SVD sur une géométrie à trois cercles intersectés. 

La décomposition de domaine de Schwarz est une méthode d’approximation d'équations aux dérivées partielles, consistant à scinder notre domaine en plusieurs sous-domaines, sur lesquels le système différentiel est approximé successivement, avec par exemple une méthode d'éléments finis de type Galerkin.

Cette méthode est en réalité assez peu utilisée, du fait de sa lenteur de convergence. Il est possible d'exploiter une propriété des erreurs d'itérations, notamment sur les bords avec la méthode d'accélération d'Aitken. Cependant, même cette méthode d’accélération n’est pas suffisante d’un point de vue performance.
En effet, il est parfois nécessaire de faire un grand nombre d'itérations de la décomposition de Schwarz pour déterminer correctement l’opérateur linéaire d’erreur. 
Une autre méthode consiste à introduire une décomposition en valeurs singulières. Cette décomposition permet de trouver les vecteurs d’erreurs les plus significatifs, et ainsi de connaître une bonne approximation de l'opérateur d'erreur au bord. C'est la méthode d'accélération d'Aitken-SVD.

Cette décomposition de Schwarz accéléré via Aitken-SVD a été montrée expérimentalement plus efficace que la décomposition de Schwarz classique. Cependant, une question subsiste : Quelles sont les performances de cette méthode quand une partie du domaine complet est incluse dans plus de deux sous-domaines ?

La méthode de Schwarz-Aitken-SVD est d'abord implémentée sur une géométrie avec deux cercles intersectés puis trois cercles intersectés, et elle est comparée à une méthode de Schwarz classique. 

Les résultats ainsi que plus de détails sont disponibles dans le rapport du projet.

Remerciements à M. Tromeur-Dervout pour ce sujet et son aide.
