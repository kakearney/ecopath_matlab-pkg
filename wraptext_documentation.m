%% |wraptext| documentation
% This function formats long strings into wrapped text of specified width. 
% 
%% Syntax
% 
%  wraptext(str) 
%  wraptext(str,width) 
%  strw = wraptext(...)
% 
%% Description 
% 
% |wraptext(str)| wraps the text in |str| and prints it in the command window. 
%
% |wraptext(str,width)| wraps text to a specified |width| in characters. If
% |width| is not specified, the current width of your command window is used.
% 
% |strw = wraptext(...)| writes output to a variable in your workspace. 
% 
%% Example 1: 
% Print text to your command window: 

str = 'Lorem interdum nescio ubi sordida iter ducis, interdum nec video cur.  Ego coniecto ego servo in ludum properamus, multa cervisia et multa vagabundos inruisse.  Facilius quam mori dum circum. Unum tempus habui amicos ma, etiam pa habebat, Et lumbare caedebat ostendam exclamavit, Nympha iubet operam me descendit vetusto Tennessium, Facilius quam mori dum circum.  In aetatem veni inuenit puellam Tuscaloosa bar, Ipsa purgata a me foras, et ledo eam clanculum, Et volebat occidere dolore emi vino tripudio agmine, Facilior visa est quam iustus expectans circa mori.  Amicus quidam dixit se scire facile pecuniae, Et vir fratrem suum fecerunt depraedantes avolavimus, Posse me adsecutus est veneno ad me Muskogee, Suus duobus annis circiter mori iustus expectans. Quod de carcere catenisque interdum Im Possedi amicum tandem, Ille ne occidas ne fureris aut fallere aut in potu aut in mendacio, Codeine nomen suum, quod ipse vidi Qui viatoribus, Simul sumus agnus dei haerere et morietur.';
wraptext(str) 

%% Example 2: 
% Write a string to your workspace and make it does not exceed 33
% characters width: 

str = 'Lorem interdum nescio ubi sordida iter ducis, interdum nec video cur.  Ego coniecto ego servo in ludum properamus, multa cervisia et multa vagabundos inruisse.  Facilius quam mori dum circum. Unum tempus habui amicos ma, etiam pa habebat, Et lumbare caedebat ostendam exclamavit, Nympha iubet operam me descendit vetusto Tennessium, Facilius quam mori dum circum.  In aetatem veni inuenit puellam Tuscaloosa bar, Ipsa purgata a me foras, et ledo eam clanculum, Et volebat occidere dolore emi vino tripudio agmine, Facilior visa est quam iustus expectans circa mori.  Amicus quidam dixit se scire facile pecuniae, Et vir fratrem suum fecerunt depraedantes avolavimus, Posse me adsecutus est veneno ad me Muskogee, Suus duobus annis circiter mori iustus expectans. Quod de carcere catenisque interdum Im Possedi amicum tandem, Ille ne occidas ne fureris aut fallere aut in potu aut in mendacio, Codeine nomen suum, quod ipse vidi Qui viatoribus, Simul sumus agnus dei haerere et morietur.';
wrapped_string = wraptext(str,33)

%% Author Info: 
% This function was written by <http://www.chadagreene.com Chad A. Greene> of the University of Texas
% at Austin's Institute for Geophysics (UTIG), September 2015. 