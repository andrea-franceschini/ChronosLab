function param = parse_xml(file_XML,ntv0)

node = xml2struct(file_XML);
if isfield(node,'ChronosLab')
   node = getfield(node,'ChronosLab');
   if isfield(node,'problem')
      problem = getfield(node,'problem');
      problem = problem.Attributes;
      if isfield(problem,'type')
         prob_type = getfield(problem,'type');
         if ( ~( strcmpi(prob_type,'mech') || strcmpi(prob_type,'lapl') ) )
            error('Wrong problem type: %s\n',prob_type);
         end
      else
         error('problem type must be defined.');
      end
   else
      error('problem must be defined.');
   end
   if isfield(node,'preconditioner')
      node = getfield(node,'preconditioner');
      node = node.Attributes;
   end
   if isfield(node,'smooth_type')
      smooth_type = getfield(node,'smooth_type');
   else
      smooth_type = 'mediumFSAI';
   end
   param.amg = parser_amg();
   param.tspace = parser_tspace(node,ntv0);
   param.smoother = parser_smoother(node,smooth_type);
   param.coarsen = parser_coarsen(node,prob_type);
   param.prolong = parser_prolong(node,prob_type);
else
   error('ChronosLab must be defined.');
end

end
