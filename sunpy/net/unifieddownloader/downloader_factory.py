from sunpy.util.datatype_factory_base import BasicRegistrationFactory
from sunpy.net.attr import *
from sunpy.net.vso.attrs import *
from sunpy.util.datatype_factory_base import NoMatchError
from sunpy.util.datatype_factory_base import MultipleMatchError
from sunpy.net.unifieddownloader.client import GenericClient
from sunpy.util import print_table

__all__ = ['UnifiedDownloader']

class UnifiedResponse(list): 

    def __init__(self, lst):
        
	tmplst = []
        for block in lst:
	    block[0].client = block[1]
	    tmplst.append(block[0])
	super(UnifiedResponse, self).__init__(tmplst)
    
    def __len__(self):
        
	ans = 0
	for qblock in self:
	    ans += len(qblock)
	return ans
    
    def __str__(self):
        
	table =[
	        [
		     qrblock.time.t1.strftime('%Y/%m/%d'),
		     qrblock.time.t2.strftime('%Y/%m/%d'),
		     qrblock.source,
		     qrblock.instrument,
		     qrblock.url
		]
		for block in self for qrblock in block
	       ]
	table.insert(0,['----------', '--------', '------', '----------', '---'])
	table.insert(0,['Start time', 'End time', 'Source', 'Instrument', 'URL'])

        return print_table(table, colsep = '  ', linesep = '\n')


qwalker = AttrWalker()

@qwalker.add_creator(AttrAnd)
def _create(wlk, query, dobj):
    #qwalker calls this function on finding AttrAnd object in query.
    qresponseobj, qclient = dobj._get_registered_widget(*query.attrs)
    return [(qresponseobj, qclient)]


@qwalker.add_creator(AttrOr)
def _create(wlk, query, dobj):    
    #qwalker calls this function on finding Attror object in query.
    qblocks = []
    for iattr in query.attrs:
        qblocks.extend(wlk.create(iattr, dobj))
    
    return qblocks


class UnifiedDownloaderFactory(BasicRegistrationFactory):


    def query(self, *query):
        '''
        and_ tranforms query into disjunctive normal form 
	ie. query is now of form A & B or ((A & B) | (C & D))
	This helps in modularising query into parts and handling each of the parts individually.
	Input:
	query: VSO style query.Attributes from JSOC,VSO both can be used.
	output: List of tuples of form(queryresponse,instance of selected client).
        '''
	query = and_(*query)	
        return UnifiedResponse(qwalker.create(query, self))

    def get(self, qr, **kwargs):
        '''
	Downloads the data.
	Input:
	List of tuples of form(queryresponse,instance of selected client).
	Output: 
	List of Results objects returned by individual clients
	'''
	reslist =[]
    	for block in qr:
		reslist.append(block.client.get(block, **kwargs))
	
	return reslist

    def __call__(self, *args, **kwargs):
        pass


    def _check_registered_widgets(self, *args, **kwargs):
        '''Factory helper function'''
        candidate_widget_types = list()
	for key in self.registry:
	    
	    if self.registry[key](*args):
	        candidate_widget_types.append(key)
            
	n_matches = len(candidate_widget_types)
	if n_matches == 0:
	    if self.default_widget_type is None:
	        raise NoMatchError("Query {0} can not be handled in its current form".format(args))
	    else:
	        return  [self.default_widget_type]
	elif n_matches > 1:
	    for candidate_client in candidate_widget_types:
	        if issubclass(candidate_client, GenericClient):
	     		return [candidate_client]
	    raise MultipleMatchError("Too many candidates clients can service your query {0}".format(args))

        return candidate_widget_types

    def _get_registered_widget(self, *args, **kwargs):
        '''Factory helper function'''
        candidate_widget_types = self._check_registered_widgets(*args)
	tmpclient = candidate_widget_types[0]()
	return tmpclient.query(*args), tmpclient


UnifiedDownloader = UnifiedDownloaderFactory(additional_validation_functions = ['_can_handle_query'])


