The breaking waves simulations are using the CFD toolbox waves2Foam. 

For waves2Foam installation, users can refer to the link: https://openfoamwiki.net/index.php/Contrib/waves2Foam.

Please not that when running breaking wave simulation cases, it can occasionaly occur that the simulation blows up at some point. 
This is due to numerical discretization problems which can be solved by reconstructing a near time step (e.g. one wave period earlier) 
and decomposing the domain from this time step. Then resume the simulation from here. 
