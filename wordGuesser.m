%**************************************************************%
% @BeginVerbatim
% Title: Genetic Algorithm-Based Word Guessing Game
% Description: A Word Guessing Game that shows Genetic Algorithm
%              in action.            
% Version: 1.0
% Author : Gensaya, Carl Jerwin F.
% Yr&Sec : BSCS 3-3
% @EndVerbatim
%**************************************************************%

%**************************************************************%
% @function : wordGuesser
% @desc     : Entry point of the program.
% @params   : none.
% @return   : none.
%***************************************************************%
function wordGuesser()
    
    clc
    format long g
    
    %*********************@Intialization*****************************%
    answer         = input('Enter word to be guessed: ', 's');
    generation     = input('Enter number of generations: ');
    ansLength      = length(answer);
    parents        = generateParents(answer); % Populates 1st Generation
    bestOffspring  = []; % Best Offspring per Gen. Holder
    bestOffFitness = []; % Best Offspring Fitness per Gen. Holder
    generationCtr  = 1:generation;
    %**************************************************************** %
    
    %*************************@Code Body*****************************%
    % Generational Loop
    for u8Ctr = 1 : generation
        
        % Populates fitness array by using gameMaster() function.
        for u8Ctr2 = 1:size(parents)
            fitness(u8Ctr2) = gameMaster(parents(u8Ctr2,:), answer);
        end
        
        % Gets and displays best fit offspring in the population.
        [minFit, minIndex] = min(fitness);
        
        % Flag for fitter offspring.
        betterFitFlag = 0; 
        
        % Loop for offspring creation.
        while  betterFitFlag == 0
            
            % Parent Selection Index Getter
            selParIndex = parentSelectionFunction(fitness);
            
            % Performs offspring production (Crossover / Mutation) 
            offspring = guesser([parents(selParIndex(1),:);...
                        parents(selParIndex(2),:)], ansLength, answer);
            
            % Checks if offspring's fitness 
            % >= than min() fitness of its parents.
            if min(gameMaster(parents(selParIndex(1),:), answer),...
                   gameMaster(parents(selParIndex(2),:), answer))...
                 > gameMaster(offspring, answer) ...
                 || gameMaster(offspring, answer) == 0;
             
               % Replaces parent with the highest fitness.
               [~, maxIndex] = max(fitness);
               parents(maxIndex,:) = offspring; 
               bestOffspring = [bestOffspring; offspring];
               bestOffFitness = [bestOffFitness; gameMaster(offspring, answer)];
               
               % Ticks flag to break while loop.
               betterFitFlag = 1; 
            end
        end
    end
    %****************************************************************%
    
    %*********************@Display_Output****************************%
    % Displays Word Guesser Stats Table
     Generation = reshape(generationCtr, [generation, 1]);
     BestGuess = cellstr(bestOffspring);
     CostValue = bestOffFitness;
     WordGuesserStat = table(Generation, BestGuess, CostValue)
     
     % Plots Cost Value Vs Generation
     figure('Name','Word Guesser using Genetic Algorithm',...
            'NumberTitle','off');
     plot(Generation, CostValue, '-o', 'LineWidth', 1.5);  
     title('Cost Vs Generation');
     xlabel('Number of Generations');
     ylabel('Cost Value');
    %****************************************************************%
end

%**************************************************************%
% @function : gameMaster
% @desc     : Computes for the fitness of a word.
% @params   : guess  - word to be evaluated 
%             answer - word to be Guessed
% @return   : fitness - computed fitness
%***************************************************************%
function fitness = gameMaster(guess, answer)

    %*********************@Intialization*****************************%
    fitness = 0;
    %****************************************************************%
    
    %*************************@Code Body*****************************%
    % Loops through every character to compute for fitness.
    for u8Ctr = 1:length(guess)
        fitness = fitness + (guess(u8Ctr) - answer(u8Ctr))^2;
    end
    %****************************************************************%
    
end

% Generates Starting Population
function parents = generateParents(answer) 

    %*********************@Intialization*****************************%
    charSet    = char(['a':'z' 'A':'Z']); % Initialize charSet.
    parents    = [];
    noOfParents = 10;
    %****************************************************************%
    
    %*************************@Code Body*****************************%
    for u8Ctr = 1 : noOfParents
        randomize  = ceil(length(charSet)*rand(1,length(answer)));
        convLetter = charSet(randomize);
        parents = [parents; convLetter];
    end
    %****************************************************************%
    
end

% Parent Selection Function
function selParIndex = parentSelectionFunction(fitness)

    %*********************@Intialization*****************************%
    [~, sortedIndex] = sort(fitness(:), 'ascend');
    %****************************************************************%
    
    %*************************@Code Body*****************************%
    selParIndex = sortedIndex(1:2);
    %****************************************************************%
    
end

%**************************************************************%
% @function : guesser
% @desc     : Performs crossover and mutation to offspring.
% @params   : parents   - selected parents
%             ansLength - length of the word to be guessed
%             answer    - word to be guessed
% @return   : offspring - offspring generated by the GA
%***************************************************************%
function offspring = guesser(parents, ansLength, answer)
    
    % Single Point Crossover
    %*********************@Intialization*****************************%
    parentA     = parents(1,:);
    parentB     = parents(2,:);
    startInd    = randi([1, ansLength - 1]);
    endInd      = randi([startInd, ansLength]);
    firstCross  = startInd:endInd;
    secondCross = setdiff(1:ansLength, firstCross);
    %****************************************************************%
    
    %*************************@Code Body*****************************%
    if gameMaster(parentA, answer) < gameMaster(parentB, answer)
       if length(firstCross) > length(secondCross)
            for u8Ctr = 1 : length(firstCross)
                crossOffspring(firstCross(u8Ctr)) =...
                    parentA(firstCross(u8Ctr));
            end
            
            for u8Ctr = 1 : length(secondCross)
                crossOffspring(secondCross(u8Ctr)) =...
                    parentB(secondCross(u8Ctr));
            end
       else
            for u8Ctr = 1 : length(firstCross)
                crossOffspring(firstCross(u8Ctr)) =...
                    parentB(firstCross(u8Ctr));
            end
            
            for u8Ctr = 1 : length(secondCross)
                crossOffspring(secondCross(u8Ctr)) =...
                    parentA(secondCross(u8Ctr));
            end    
       end
       
    else
        if length(firstCross) > length(secondCross)
            for u8Ctr = 1 : length(firstCross)
                crossOffspring(firstCross(u8Ctr)) =...
                    parentB(firstCross(u8Ctr));
            end
            
            for u8Ctr = 1 : length(secondCross)
                crossOffspring(secondCross(u8Ctr)) =...
                    parentA(secondCross(u8Ctr));
            end
       else
            for u8Ctr = 1 : length(firstCross)
                crossOffspring(firstCross(u8Ctr)) =...
                    parentA(firstCross(u8Ctr));
            end
            
            for u8Ctr = 1 : length(secondCross)
                crossOffspring(secondCross(u8Ctr)) =...
                    parentB(secondCross(u8Ctr));
            end    
       end
    end
    %****************************************************************%
    
    % Mutation
    %*********************@Intialization*****************************%
    
    % Initialize Characted Set to be Used
    charSet    = char(['a':'z' 'A':'Z']); 
    randLetter = randi([1, 52]);
    randIndex  = randi([1, length(crossOffspring)]);
    convLetter = charSet(randLetter);
    
    %****************************************************************%
    
    %*************************@Code Body*****************************%
    crossOffspring(randIndex) = convLetter;
    offspring = crossOffspring;
    %***************************************************************%
    
end