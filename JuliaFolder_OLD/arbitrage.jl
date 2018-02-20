# This is a script to determine whether there is an arbitrage opportunity for a
# set of bookmaker's odds
# In particular this is for a game which can either be won or lost eg tennis

num_corporates = 5

odds = zeros(2,num_corporates)
true_odds = [1.558; 2.460]
true_percentage = 1/true_odds[1] + 1/true_odds[2]

for i=1:num_corporates
  perturbation_odds_1 = (rand()-0.5)/(true_odds[1]*3)
  #perturbation_odds_2 = odds[2,i] - (odds[1,i]+perturbation_odds_1)/(true_percentage*odds[1] + true_percentage*perturbation_odds_1 - 1)
  odds[1,i] = true_odds[1] + perturbation_odds_1
  odds[2,i] = true_odds[2] - perturbation_odds_1
end

percentages = zeros(1,num_corporates)
for i=1:num_corporates
  percentages[i] = 1/odds[1,i] + 1/odds[2,i]
end

best_odds = zeros(2,1)
best_odds[1] = maximum(odds[1,:])
best_odds[2] = maximum(odds[2,:])

percentages_best_odds = 1/best_odds[1] + 1/best_odds[2]

println(best_odds)
println(percentages_best_odds)

x = 1
#y = 0.5*(x*best_odds[1] - x + x/(best_odds[2]-1))
#y = (best_odds[1]+best_odds[2])/(1+best_odds[2])
y = best_odds[1]*x/best_odds[2]
profit_if_1 = best_odds[1]*x - x - y
profit_if_2 = best_odds[2]*y - y - x

expected_profit = profit_if_1/(percentages_best_odds*best_odds[1]) + profit_if_2/(percentages_best_odds*best_odds[2])
percetage_return = expected_profit/(x+y)
println(percetage_return)
