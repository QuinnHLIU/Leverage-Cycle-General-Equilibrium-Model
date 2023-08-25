# Leverage-Cycle-General-Equilibrium-Model
example code for leverage cycle models

The project firstly consists several slides briefly introducing the collateral equilibrium model by Fostel and Geanakoplos (2008, 2012). The code folder consists sample matlab function to showcase the basic idea.

## Section 1. The Leverage Economy and Concept of Collateral Equilibrium

The model considers a simple two period model with a safe asset and a risky asset. Agents have heterogenous belief about future asset prices, and the equilibrium should be determined by (i) a marginal buyer of the risky asset (with leverage) (ii) the equilibrium asset prices.

The equilibrium is solved by equating both the indifference condition of the marginal buyer, and the market clearing condition.

\\

Code: **Leverage.m** 

It firstly computes a no borrowing benchmark, and compare the result with "Leverage Economy" allowing for borrowing. Debts are all non-recourse hence the No-default theorem applies (Geanakoplos, 1997), i.e., only the min-max contracts will be traded.


## Section 2. The Anxious Economy

Next we consider an "anxious economy" in which there are three periods, hence introduce dynamics to asset prices. The model shows (i) debts are endogenously short-term (ii) leverage amplifies asset price fluctuations (iii) bad news on one specific asset can have contagious effect in the sense that the prices of other assets go down as well even without any news about their fundamentals.

\\

Code: **Anxious.m**


## Section 3. Financial Innovation in Multi-States

We consider financial innovation when the economy has three states. An economy with N states generally allows for N-1 types of contracts to be traded (Phelan and Toda, 2019) in the equilibrium (the min-max principle still applies, yet of course not No-default anymore). Three types of financial innovation are considered.
\\

(1) Tranching: creating a down tranche is to split the asset payoffs to an Arrow Up and an Arrow Down. This deviates from the collateral equilibrium in an important direction: allows for state-contingent promises.

(2) CDS: essentially, a tranching of Cash (Fostel and Geanakoplos, 2012a).

(3) Securitization (Pyramiding): distinct from Tranching and CDS, which is creating new promises based on existing collateral, securtization is creating new collateral --- allowing for promises to be backed by existing _promises_ with certain extent of contingency.

The more number of states, the richer the equilibrium contract space could do. The wonderful extension to continuous state space is done by Simsek (2013).
> It highlights an important mindset that, in an economy with more than two states, the "nature" of beliefs is crucial to the analysis: are you optimistic about upside to happen, or are you optimistic about downside is less likely?


\\

Code: **innovation.m**


## Section 4. Global Collateral

To consider the model in an open-economy setting, we assume two economies (Foreign and Home) with different degree of financial development. Fostel, Geanakoplos, and Phelan (2019): Home (the US) is in a Tranching Economy, while Foreign (the EU) is in a Leverage Economy.

Then we consider if the asset markets are freely opened, so both the Home and Foreign assets can be traded by Home and Foreign investors. The analysis shows that, free trade cannot eliminate the asset price dispersion, since the _collateralizability_ of these two assets still differs. As a matter of fact, international capital flows enlarges this dispersion because the Home assets are more desirable for all investors. The Foreign assets become even less valuable if Foreign investors can make tranches instead of just promising state-contingent debts.

\\

Code: **globalcollateral.m**


**References**: 
Fostel and Geanakoplos (2008): ``Leverage Cycles and the Anxious Economy," _American Economic Review_.

Fostel and Geanakoplos (2012a): ``Why Does Bad News Increase Volatility and Decrease Leverage?" _Journal of Economic Theory_.

Fostel and Geanakoplos (2012b): ``Traching, CDS, and Asset prices: How Financial Innovation Can Cause Bubbles and Crashes," _American Economic Journal of Economic Theory_.

Fostel, Geanakoplos, and Phelan (2019): ``Global Collateral and Capital Flows," NBER Working Paper.

Phelan and Toda (2019): ``Securitized markets, international capital flows, and global welfare,", _Journal of Financial Economics_.

Simsek (2013): ``Belief Disagreements and Collateral Constraints,", _Econometrica_.
