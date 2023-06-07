## Primo incontro
- Perchè le forze Fx e Fy sull'asse vengono calcolate senza moltiplicare tutte le forze per la matrice di rotazione? (R159 Simulink-Vehicle Model)
- Perchè la Delta Fx viene negativa (anche se i valori sono piccoli) quando stanno facendo una curva a sinistra e sto decellerando?
- Come mantenere un'accelerazione laterale costante?
- Nostra idea per fare il grafico dell'axle characteristic: simulare una serie di situazioni di sterzo/pedale. Da ogni simulazione estrarre le coppie slip/forza corrispondenti ad un certo valore di accelerazione. Plottare su un grafico le curve slip/forza per i diversi livelli di accelerazione (eg con colori diversi). Il problema potrebbe essere come scegliere i test case da svolgere, in modo da trovare una serie di dadti tra loro sensati.
- Perche se sterzo verso sinistra ho una velocità laterale negativa?
- Perchè la simulazione non funziona mettendo i nostri coefficienti?
- Il toe deve essere inserito solo come offset nel perfect akermann? Dobbiamo tenerne anche in considerazione nel calcolo di alpha_f e alpha_r?
- Per il camber dobbiamo modificare anche il modello, oppure è già considerato dentro la funzione camberModel? 


- Change the sign of all the alpha
- Positive camber increases the paek of the lateral force. Signo fo the camber according to the reference frame used in the wheel (pos sx, neg dx)
- AXLE CHARACTERISTIC: Dare input di pedale, sterzare con un certo angolo, la macchina comincerà a girare in torno. Una volta raggiunta una condizione stazionaria, estrarre lo slip corrispondente e i valori delle forze (che verranno poi plottati su un grafico). Ripetere questo test per diversi valori di pedale e sterzo in modo da ottenere una serie di punti sul grafico
- HANDLING DIAGRAM: curvare poco, tenedo una velocità di rotazione costante (eventualmente inerendo un controllore per la velocità longitudinale), plottare poi la velocità di avanzamento contro l'handling
- Simulazione per camber costante già mimplementata, possibilità di aggiungere la dinamica della sospensione. In particolare la variazione del camber è proporzionale all'angolo di roll. 
- Attenzione alla convenzione dei segni per l'angolo di toe
- Cambiare i segni dei self aligning torque