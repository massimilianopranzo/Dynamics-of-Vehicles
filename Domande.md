## Primo incontro
- Perchè le forze Fx e Fy sull'asse vengono calcolate senza moltiplicare tutte le forze per la matrice di rotazione? (R159 Simulink-Vehicle Model)
- Perchè la Delta Fx viene negativa (anche se i valori sono piccoli) quando stanno facendo una curva a sinistra e sto decellerando?
- Come mantenere un'accelerazione laterale costante?
- Nostra idea per fare il grafico dell'axle characteristic: simulare una serie di situazioni di sterzo/pedale. Da ogni simulazione estrarre le coppie slip/forza corrispondenti ad un certo valore di accelerazione. Plottare su un grafico le curve slip/forza per i diversi livelli di accelerazione (eg con colori diversi). Il problema potrebbe essere come scegliere i test case da svolgere, in modo da trovare una serie di dadti tra loro sensati.
- Perche se sterzo verso sinistra ho una velocità laterale negativa?
- Perchè la simulazione non funziona mettendo i nostri coefficienti?
- Il toe deve essere inserito solo come offset nel perfect akermann? Dobbiamo tenerne anche in considerazione nel calcolo di alpha_f e alpha_r?
- Per il camber dobbiamo modificare anche il modello, oppure è già considerato dentro la funzione camberModel? 
