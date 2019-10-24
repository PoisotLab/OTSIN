# Modeling mutualistic networks using optimal transport

Mutualistic networks, such as plant-pollinator, bird-seed dispersal, or plant-ant networks, can be understood as systems that try to maximize both utility (every species tends to its favorite partners) and freedom of choice in interactions. I will derive this from first principles.[^3f30df83]

[^3f30df83]: I think this should not necessarily be limited to mutualistic networks. Large part of the theory should also be relevant when only one of the types of species tries to maximize utility, e.g. host-parasite and food webs.

1. There are two types of ecological partners, e.g. pollinator species and plant species. Let us call the types $A$ (e.g. animal) and $B$ (e.g. plant) and there are $m$ and $n$ species of each type respectively.
2. There is some (hypothetical) resource or *currency* exchanged by those partners, in the case of pollination, this might be nectar or pollen. Let the relative uptake for each species of the animals be a normalized histogram $\mathbf{a}$, i.e. $\mathbf{a}\in\mathbb{R}^m$ satisfying $a_i\geq 0$ and $\sum_ia_i=1$. Likewise, the individual plant species provide this resource distributed according to the histogram $\mathbf{b}\in\mathbb{R}^n$ satisfying $b_j\geq 0$ and $\sum_jb_j=1$.
   - If this currency is not known, it can be estimated by the relative number of visits: $$a_i = \frac{\text{number of visits performed by animal }i}{\text{total number of visits performed by all animals}}$$ and visa versa for the plant species. These values can be improved by weighting the visits with the relative uptake or production of resources by the species.
3. There exists a *coupling* between the two types of species, $Q$, an $m \times n$ matrix describing the fraction of currency each animal takes from each plant. Permissible couplings should be in agreement with the histograms of currency uptake and production. These couplings should be an element of the transport polytope of $\mathbf{a}$ and $\mathbf{b}$: $$U(\mathbf{a}, \mathbf{b}) = \{Q\in\mathbb{R}_+^{m\times n}\mid \sum_jQ_{ij}=a_i, \sum_iQ_{ij}=b_j\}\,.$$ In other words, only positive real matrices for which the row sums and columns sums are $\mathbf{a}$ and $\mathbf{b}$, respectively.
   - We can estimate the coupling as proportional with the number of visits: $$Q_{ij} = \frac{\text{number of visits performed by animal $i$ to plant $j$}}{\text{total number of visits performed by all species}}\,.$$ Again, the visits can be weighted by relative resource production/uptake.
4. Interaction between some animals and plants are more efficient or stronger preferred than others. Let $M$ be an $m \times n$ matrix representing the *utility* of the different interactions. For example, element $M_{ij}$ is the utility between animal $i$ and plant $j$. The average utitly of a system is given by $\langle M, Q\rangle_F = \sum_{ij} M_{ij}Q_{ij}$. **Ecosystems tends to maximize the average utility in short term by generating an optimal coupling, i.e. species choose their interactions to increase global utility.**
   - Without loss of generality, let us assume that there utilities are all values between 0 and 1. This can always be obtained by scaling and translating arbitrary units of currency.
5. Couplings are also driven by random or *neutral* processes. This can be quantified by the *entropy* of a coupling: $$H(Q) = -\sum_{i,j}Q_{ij}\log Q_{ij}\,.$$ **Ecosystems tend to increase the entropy of the couplings by random processes and incomplete information**.
   - The coupling that maximizes the entropy is computed as $\mathbf{a}\mathbf{b}^\top$.
6. There is a trade-off between maximizing the average utility of an ecosystem and the entropy of the couplings. This is quantified by a positive parameter $\lambda$. As such, the coupling of an ecosystem for given animal and plant distributions is obtained by solving the following optimization problem: $$ \max_{Q\in U(\mathbf{a}, \mathbf{b})} \langle M, Q\rangle_F+ \lambda H(Q)\,.$$
   - This is a *strictly convex* optimization problem. There is always an unique solution!
   - The optimal coupling $Q^\star$ can be found quickly using *Sinkhorn's algorithm* [@Cuturi2013].
   - If $\lambda \rightarrow 0$, the system is only driven by maximizing the utility of species. The coupling can be found reasonable efficiently using dynamic programming.
   - If $\lambda \rightarrow \infty$, the system only driven by neutral processes. In this case $Q^\star = \mathbf{a} \mathbf{b}^\top$.
7. At short time scales, the distribution of the species remain fixed with only the interactions changing to maximize utility and entropy. **At longer time scales, we assume that the distributions of the species can also change to increase both utility and random interactions.** Hence, on the long run, the optimal coupling between the species is obtained by solving: $$\max_{Q \text{ s.t. } Q_{ij}\geq 0, \sum_{i,j}Q_{ij}=1} \langle M, Q\rangle_F+ \lambda H(Q)\,.$$
   - The optimal coupling on the long term, when the system is in equilibrium is given by $Q_{ij}^\triangle=\frac{\exp{M_{ij}/\lambda}}{\sum_{kl}\exp{M_{kl}/\lambda}}$. An appropriate scaling factor needs to be chosen such that $Q^\triangle$ is normalized.
      - This is the *softmax* of the matrix $M/\lambda$, there is a long history of this in decision theory, reinforcement learning and game theory.
   - Given the observed coupling matrix $P$, we can fit the utility by minimizing the *Kullback-Leibler* divergence $$D_\text{KL} = \sum_{ij} P_{ij}\log(\frac{Q_{ij}}{P_{ij}})\,.$$
      - If certain interactions do not occur, this will result in an $M_{ij}→-\infty$. This is consistent with what we expect.
      - In practice, $P$ is only an estimate, we can add *pseudocounts* or regularization to $M$ to correct for potential false negatives in our network.
      - $M$ is invariant for symmetries such as scaling or translation. Linear constraints on $M$ can be added.
8. Given that one of the species distributions is fixed, we can also easily obtain the optimal distribution maximizing utility (with a minimal entropy). Suppose that we want to sow different plants optimal with respect to pollinator community, this is the following optimization problem: $$\max_{Q \text{ s.t. } Q_{ij}\geq 0, \sum_{j}Q_{ij}=a_i} \langle M, Q\rangle_F+ \lambda H(Q)\,,$$ which has an optimal solution in the form of $$Q^\star(\mathbf{a})=\beta_ie^{M_{ij}/\lambda}\,,$$ where the $\beta_j$ parameters are chosen such that $Q^\star$ matches the distribution of the given pollinators.

## Summary & thoughts

- Using optimal transport, we can simulate the couplings between animals and plants in the short term.
- Adding a constant to the utility: $M'=M+c$ (i.e. setting the baseline utility) has no effect on the optimal coupling. We set the lowest value in the matrix to zero by definition.
- Likewise, scaling the utility $M'=cM$ also has no effect on the optimal couplings if the entropy parameter is scaled proportionally.
- Assuming that the species distributions are in an equilibrium, i.e. species abundances are such to maximize both average utility and entropy, we can estimate the utilities and entropy parameter from an co-occurrence matrix, up to a linear rescaling.
- systems with small $\lambda$ are more *efficient* while systems with large $\lambda$ are more *stable*.
- There is a strong link with statistical physics.
- Link this to a dynamical model? E.g. replicator equation?
- See simulation experiments on different types of network structures.


## Literature survey

- **Optimal foraging theory** is a behavioral ecology model that helps predict how an animal behaves when searching for food. Although obtaining food provides the animal with energy, searching for and capturing the food require both energy and time. To maximize fitness, an animal adopts a foraging strategy that provides the most benefit (energy) for the lowest cost, maximizing the net energy gained. OFT helps predict the best strategy that an animal can use to achieve this goal.
- MaxEnt principle [@Harte2014]
- @Magrach2017: Honeybee spillover article.
- [The case for neutral ecology](https://www.sciencedirect.com/science/article/pii/S0169534712000237?dgcid=raven_sd_recommender_email)
- @Urban-Mead2017 performed an experiment on pollinator visits, removed certain types of flowers (open vs tubular) in a plots (control vs treatment) and observed how the network changed.
- @Brosi2017 measured interaction networks (number of visits) with and without the exclusion of a type of bumblebee.
- Pollination networks from natural and anthropogenic-novel communities show high structural similarity
- @Albrecht2018: nice plant-frugivore dataset with traits and phylogeny. Might illustrate a regression model?
- @Santamaria2007: *Linkage rules for plant-pollinator networks: Trait complementarity or exploitation barriers?*
    - forbidden links/trait complementarity vs neutral hypothesis
        - forbidden links split in traits and phylogeny
        - mix of both works best
    - 'deterring floral parasites is at least as important as increasing pollination efficiency in the evolution of plant–pollinator networks'
    - interessante discussie neutral hypothesis
- @Cirtwill2019: The authors propose a Bayesian model to provide a posteriory on distribution for species interactions. It is mainly useful when there are several networks, where some species do not meet. Uncertainty is divided into interaction probability, co-occurence proability and detectabitlity and the model takes all these things into account.
