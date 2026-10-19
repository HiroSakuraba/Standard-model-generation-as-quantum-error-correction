-- Move \tag{n} out of display math in HTML output so equation numbers survive MathML.
function Math(el)
  if el.mathtype ~= 'DisplayMath' then return nil end
  local n = el.text:match('\\tag{([^}]*)}')
  if not n or not FORMAT:match('html') then return nil end
  local stripped = el.text:gsub('\\tag{[^}]*}', '')
  return { pandoc.Span({ pandoc.Str('(' .. n .. ')') }, { class = 'eqno' }), pandoc.Math('DisplayMath', stripped) }
end
